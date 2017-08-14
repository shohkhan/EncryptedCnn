#include <vector>
#include "/HElib/HElib-master/src/FHE.h"
#include "/HElib/HElib-master/src/timing.h"
#include "/HElib/HElib-master/src/EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>
#include <fstream>
#include <sstream>
#include <sys/time.h>
#include <omp.h>
#include <pthread.h>
#include <thread>
#include <mutex>
using namespace std;
mutex mu;

vector<int> filterPlainImage(int stridex, int stridey, int filterSize,
		int width, int length, vector<int> plainImage, vector<int> filter) {
	int pad = filterSize / 2;
	vector<int> filteredPlainImage;

	for (int l = pad; l < length - pad; l += stridey) {
		for (int e = pad; e < width - pad; e += stridex) {
			int v = 0;
			for (int i = 0; i < filterSize; i++) {
				int y = l - pad + i;
				if (y >= 0 && y < length) {
					for (int j = 0; j < filterSize; j++) {
						int x = e - pad + j;
						if (x >= 0 && x < width) {
							int v1 = plainImage[y * width + x];
							int mul = filter[i * filterSize + j];
							v1 *= mul;
							v += v1;
						}
					}
				}
			}
			filteredPlainImage.push_back(v);
		}
	}
	return filteredPlainImage;
}

void filterEncryptedImage(int stridex, int stridey, int filterSize, int width,
		int length, vector<Ctxt> encImage, vector<int> filter, Ctxt zero,
		vector<vector<Ctxt>> & filteredImages, int sequence) {
	//this_thread::sleep_for(chrono::seconds(2));
	int pad = filterSize / 2;
	vector<Ctxt> filteredImage;

	for (int l = pad; l < length - pad; l += stridey) {
		for (int e = pad; e < width - pad; e += stridex) {
			Ctxt v = zero;
			for (int i = 0; i < filterSize; i++) {

				int y = l - pad + i;
				if (y >= 0 && y < length) {
					for (int j = 0; j < filterSize; j++) {
						int x = e - pad + j;
						if (x >= 0 && x < width) {
							Ctxt v1 = encImage[y * width + x];
							int mul = filter[i * filterSize + j];
							if (mul != 0) {
								v1.multByConstant(to_ZZ(mul));
								v.addCtxt(v1);
								v1.clear();
							}
						}
					}
				}
			}
			filteredImage.push_back(v);
			v.clear();
		}
	}
	//mu.lock();
	filteredImages[sequence] = filteredImage;
	//mu.unlock();
	filteredImage.clear();
}

void filterEncryptedImageAdd(int stridex, int stridey, int filterSize, int width,
		int length, vector<Ctxt> encImage, vector<int> filter, Ctxt zero,
		vector<vector<Ctxt>> & filteredImages, int sequence) {
	int pad = filterSize / 2;
	int a = 0;
	for (int l = pad; l < length - pad; l += stridey) {
		for (int e = pad; e < width - pad; e += stridex) {
			Ctxt v = zero;
			for (int i = 0; i < filterSize; i++) {

				int y = l - pad + i;
				if (y >= 0 && y < length) {
					for (int j = 0; j < filterSize; j++) {
						int x = e - pad + j;
						if (x >= 0 && x < width) {
							Ctxt v1 = encImage[y * width + x];
							int mul = filter[i * filterSize + j];
							if (mul != 0) {
								v1.multByConstant(to_ZZ(mul));
								v.addCtxt(v1);
								v1.clear();
							}
						}
					}
				}
			}
			filteredImages[sequence][a].addCtxt(v);
			a ++;
			//filteredImage.push_back(v);
			v.clear();
		}
	}
	//mu.lock();
	//filteredImages[sequence] = filteredImage;
	//mu.unlock();
	//filteredImage.clear();
}

vector<int> poolPlainImage(int poolSize, int width, int length,
		vector<int> plainImage) {
	vector<int> pooledPlainImage;

	for (int l = 0; l < length; l += 2) {
		for (int e = 0; e < width; e += 2) {
			int v = 0;
			for (int i = 0; i < poolSize; i++) {
				int y = l + i;
				if (y >= 0 && y < length) {
					for (int j = 0; j < poolSize; j++) {
						int x = e + j;
						if (x >= 0 && x < width) {
							v += plainImage[(l + i) * width + e + j];
						}
					}
				}
			}
			pooledPlainImage.push_back(v);
		}
	}
	return pooledPlainImage;
}

void poolEncryptedImage(int poolSize, int width, int length, vector<Ctxt> encImage,
		Ctxt zero, vector<vector<Ctxt>> & pooledImages, int sequence) {
	vector<Ctxt> pooledImage;
	for (int l = 0; l < length; l += 2) {
		for (int e = 0; e < width; e += 2) {
			Ctxt v = zero;
			for (int i = 0; i < poolSize; i++) {
				int y = l + i;
				if (y >= 0 && y < length) {
					for (int j = 0; j < poolSize; j++) {
						int x = e + j;
						if (x >= 0 && x < width) {
							v.addCtxt(encImage[(l + i) * width + e + j]);
						}
					}
				}
			}
			pooledImage.push_back(v);
			v.clear();
		}
	}
	pooledImages[sequence] = pooledImage;
	pooledImage.clear();
}

void convOperationParallel(int stridex, int stridey, int filterSize, int width,
		int length, vector<Ctxt> encImage, int noOfFilters,
		vector<vector<int> > filters, Ctxt zero,
		vector<vector<Ctxt>> & filteredImages) {

	for (int i = 0; i < noOfFilters; i++) {
		vector<Ctxt> fi;
		filteredImages.push_back(fi);
	}
	thread *tt = new thread[noOfFilters];
	for (int i = 0; i < noOfFilters; i++) {
		tt[i] = thread(filterEncryptedImage, stridex, stridey, filterSize,
				width, length, encImage, filters[i], zero, ref(filteredImages),
				i);
	}
	for (int i = 0; i < noOfFilters; i++) {
		tt[i].join();
	}
}

void convOperationParallelAdd(int stridex, int stridey, int filterSize, int width,
		int length, vector<Ctxt> encImage, int noOfFilters, vector<vector<int>> filters,
		Ctxt zero, vector<vector<Ctxt>> & filteredImages) {
	thread *tt = new thread[noOfFilters];
	for (int i = 0; i < noOfFilters; i++) {
		tt[i] = thread(filterEncryptedImageAdd, stridex, stridey, filterSize,
				width, length, encImage, filters[i], zero, ref(filteredImages),
				i);
	}
	for (int i = 0; i < noOfFilters; i++) {
		tt[i].join();
	}
}

void convOperationParallelWithMainThread(int stridex, int stridey,
		int filterSize, int width, int length, vector<Ctxt> encImage,
		int noOfFilters, vector<vector<int> > filters, Ctxt zero,
		vector<vector<Ctxt>> & filteredImages) {

	for (int i = 0; i < noOfFilters; i++) {
		vector<Ctxt> fi;
		filteredImages.push_back(fi);
	}
	thread *tt = new thread[noOfFilters - 1];
	for (int i = 0; i < noOfFilters - 1; i++) {
		tt[i] = thread(filterEncryptedImage, stridex, stridey, filterSize,
				width, length, encImage, filters[i], zero, ref(filteredImages),
				i);
	}
	filterEncryptedImage(stridex, stridey, filterSize, width, length, encImage,
			filters[noOfFilters - 1], zero, ref(filteredImages),
			noOfFilters - 1);
	for (int i = 0; i < noOfFilters - 1; i++) {
		tt[i].join();
	}
}

void convOperationSerial(int stridex, int stridey, int filterSize, int width,
		int length, vector<Ctxt> encImage, int noOfFilters,
		vector<vector<int> > filters, Ctxt zero,
		vector<vector<Ctxt>> & filteredImages) {

	for (int i = 0; i < noOfFilters; i++) {
		vector<Ctxt> fi;
		filteredImages.push_back(fi);
	}
	for (int i = 0; i < noOfFilters; i++) {
		filterEncryptedImage(stridex, stridey, filterSize, width, length,
				encImage, filters[i], zero, ref(filteredImages), i);
	}
}

void poolOperationParallel(int poolSize, int width, int length, Ctxt zero,
		vector<vector<Ctxt>> encImages, vector<vector<Ctxt>> & pooledImages) {

	for (unsigned int i = 0; i < encImages.size(); i++) {
		vector<Ctxt> pooledImage;
		pooledImages.push_back(pooledImage);
	}
	thread *tt = new thread[encImages.size()];
	for (unsigned int i = 0; i < encImages.size(); i++) {
		tt[i] = thread(poolEncryptedImage, poolSize, width, length,
				encImages[i], zero, ref(pooledImages), i);
	}
	for (unsigned int i = 0; i < encImages.size(); i++) {
		tt[i].join();
	}
}

int main(void) {
	//For stop watch
	double times;
	const FHEtimer *tmp;
	resetAllTimers();
	setTimersOn();

	//Fixed
	//int batchSize = 2000;
	int filterSize = 3;
	int poolSize = 2;
	int stridex = 1;
	int stridey = 1;
	int maxFilterSize = 50;

	//Not fixed
	int width = 28;
	int length = 28;

	/*** BEGIN INITIALIZATION ***/
	cout << "Initiating scheme ... " << flush;
	long m = 0;		// Specific modulus
	//long p = 19777;
	long p = 19777;
	// Plaintext base [default=2], should be a prime number
	long r = 1;		// Lifting [default=1]
	long L = 6;		// Number of levels in the modulus chain [default=heuristic]
	long c = 2;		// Number of columns in key-switching matrix [default=2]
	long w = 64;	// Hamming weight of secret key
	long d = 1;		// Degree of the field extension [default=1]
	long k = 80;	// Security parameter [default=80]
	long s = 0;		// Minimum number of slots [default=0]

	long seed = 0;
	SetSeed(ZZ(seed));

	m = FindM(k, L, c, p, d, s, 0, false);// Find a value for m given the specified values
	FHEcontext context(m, p, r);			// Initialize context
	buildModChain(context, L, c);// Modify the context, adding primes to the modulus chain

	FHESecKey secretKey(context);
	const FHEPubKey& publicKey = secretKey;
	secretKey.GenSecKey(w);					// A Hamming-weight-w secret key
	ZZX G;
	G = makeIrredPoly(p, d);
	addSome1DMatrices(secretKey);// compute key-switching matrices that we need
	cout << "done" << endl;
	/*** END INITILIZATION ***/

	Ctxt zero = Ctxt(publicKey);
	cout << "# Image size: " << width * length << endl;
	/*** BEGIN: Image Encryption ***/
	vector<Ctxt> encImage;
	vector<int> plainImage;
	cout << "Encrypting image ... " << flush;
	FHE_NTIMER_START(enc);
	int singleImage = 1;
	if (singleImage == 1) {

		//cout << endl;
		for (int i = 0; i < length * width; i++) {
			Ctxt ctx1(publicKey);
			publicKey.Encrypt(ctx1, to_ZZX(i));
			encImage.push_back(ctx1);
			plainImage.push_back(i);
			//cout << i << " " << flush;
			//if ((i + 1) % 5 == 0)
				//cout << endl;
		}
		//cout << endl;
	} else {
		//Batch
		EncryptedArray ea(context, G);
		long nslots = ea.size();
		cout << "(Slots: " << nslots << ")" << flush;

		vector<long> lp;
		for (int b = 0; b < nslots; b++) {
			lp.push_back(1);
		}

		Ctxt ctx1(publicKey);
		ea.encrypt(ctx1, publicKey, lp);

		for (int i = 0; i < length * width; i++) {
			vector<long> lp1;
			for (int b = 0; b < nslots; b++) {
				lp1.push_back(1);
			}
			Ctxt c1(publicKey);
			ea.encrypt(c1, publicKey, lp1);
			encImage.push_back(ctx1);
			plainImage.push_back(i % width + 1);
		}
		ctx1.clear();

	}
	cout << "done" << endl;
	/*** END: Image Encryption ***/
	FHE_NTIMER_STOP(enc);
	tmp = getTimerByName("enc");
	times = tmp->getTime();
	cout << "Time for enc = " << times << endl;

	//ofstream myfile;
	//myfile.open("example.txt");
	//myfile << encImage;
	//myfile.close();

	/*** END: Image Encryption ***/

	cout << "generating filters ... " << flush;
	vector<vector<int> > filters;
	for (int f = 0; f < maxFilterSize; f++) {
		vector<int> filter;
		for (int i = 0; i < filterSize * filterSize; i++) {
			filter.push_back(1);
		}
		filters.push_back(filter);
	}
	cout << "done" << endl;

	/**** compare plain and encrypted ****/
	//Part 1: Pooling

	/*vector<int> pooledPlainImage;
	pooledPlainImage = poolPlainImage(poolSize, width, length, plainImage);
	vector<vector<Ctxt>> pooledImages;
	vector<Ctxt> pooledImage;
	pooledImages.push_back(pooledImage);
	FHE_NTIMER_START(po1);
	poolEncryptedImage(poolSize, width, length, encImage, zero, ref(pooledImages), 0);
	//vector<Ctxt> pooledImage = poolEncryptedImage(poolSize, width, length, encImage, zero);
	FHE_NTIMER_STOP(po1);
	tmp = getTimerByName("po1");
	times = tmp->getTime();
	cout << "Time for 1 pool = " << times << endl;
	pooledImage = pooledImages[0];
	cout << "Compare Pooled Images" << endl;
	for (unsigned int i = 0; i < pooledPlainImage.size(); i++) {
		cout << pooledPlainImage[i] << " " << flush;
	}
	cout << endl;
	cout << pooledImage.size() << endl;
	for (unsigned int i = 0; i < pooledImage.size(); i++) {
		ZZX pt;
		secretKey.Decrypt(pt, pooledImage[i]);
		cout << pt[0] << " " << flush;
	}
	cout << endl;
	pooledPlainImage.clear();
	pooledImage.clear();
	pooledImages.clear();*/

	//Part 2: Filtering
	 /*vector<int> filteredPlainImage;
	 filteredPlainImage = filterPlainImage(stridex, stridey, filterSize, width,
	 length, plainImage, filters[1]);
	 FHE_NTIMER_START(fil1);
	 vector<vector<Ctxt>> filteredImages1;
	 vector<Ctxt> filteredImage;
	 filteredImages1.push_back(filteredImage);
	 filterEncryptedImage(stridex, stridey, filterSize, width, length, encImage,
	 filters[1], zero, ref(filteredImages1), 0);
	 FHE_NTIMER_STOP(fil1);
	 tmp = getTimerByName("fil1");
	 times = tmp->getTime();
	 cout << "Time for 1 filter = " << times << endl;

//	 filteredImage = filteredImages1[0];
//	 cout << "Compare Filtered Images" << endl;
//	 for (unsigned int i = 0; i < filteredImage.size(); i++) {
//		 ZZX pt;
//	 	 secretKey.Decrypt(pt, filteredImage[i]);
//	 	 cout << pt[0] << " " << filteredPlainImage[i] << endl;
//	 }
//	 filteredImage.clear();
	 filteredPlainImage.clear();
	 filteredImages1.clear();*/
	/**** end ****/

	int noOfFilters = 20;

	cout << "Conv start ... " << flush;
	cout << "Filters: " << noOfFilters << "..." << flush;
	vector<vector<Ctxt>> convOut1;

	FHE_NTIMER_START(conv1);

	convOperationParallel(stridex, stridey, filterSize, width, length, encImage,
			noOfFilters, filters, zero, ref(convOut1));
	encImage.clear();
	FHE_NTIMER_STOP(conv1);
	tmp = getTimerByName("conv1");
	times = tmp->getTime();
	cout << "end" << endl;
	cout << "# Time for conv: " << times << endl;
	cout << "# Size after conv: " << convOut1[0].size() << endl;
	cout << "# Depth after conv: " << convOut1.size() << endl;

//	for (int i=0; i < noOfFilters; i++){
//		vector<Ctxt> filteredImage = convOut1[i];
//		for (unsigned int i = 0; i < 5; i++) {
//			ZZX pt;
//			secretKey.Decrypt(pt, filteredImage[i]);
//			cout << pt[0] << " " << flush;
//		}
//		cout << endl;
//	}

	cout << "Pool start ... " << flush;

	width = 26;
	length = 26;

	FHE_NTIMER_START(pool1);

	vector<vector<Ctxt>> poolOut1;
	//poolOperationParallel(poolSize, width, length,zero, convOut1, ref(poolOut1));

	for (unsigned int i = 0; i < convOut1.size(); i++) {
		vector<Ctxt> pooledImage;
		poolOut1.push_back(pooledImage);
	}
	thread *tt = new thread[convOut1.size()];
	for (unsigned int i = 0; i < convOut1.size(); i++) {
		tt[i] = thread(poolEncryptedImage, poolSize, width, length, convOut1[i], zero,
				ref(poolOut1), i);
		//poolEncryptedImage(poolSize, width, length, convOut1[i], ref(poolOut1), i);
	}
	for (unsigned int i = 0; i < convOut1.size(); i++) {
		tt[i].join();
	}
	convOut1.clear();
	FHE_NTIMER_STOP(pool1);
	tmp = getTimerByName("pool1");
	times = tmp->getTime();
	cout << "end" << endl;
	cout << "# Time for Pool: " << times << endl;
	cout << "# Size after pool: " << poolOut1[0].size() << endl;
	cout << "# Depth after pool: " << poolOut1.size() << endl;

//	for (int i=0; i < 5; i++){
//		vector<Ctxt> filteredImage = poolOut1[i];
//		for (unsigned int i = 0; i < 5; i++) {
//			ZZX pt;
//			secretKey.Decrypt(pt, filteredImage[i]);
//			cout << pt[0] << " " << flush;
//		}
//		cout << endl;
//	}

	noOfFilters = 50;
	cout << "Conv start ... " << flush;
	cout << "Filters: " << noOfFilters << "..." << flush;
	width = 13;
	length = 13;

	FHE_NTIMER_START(conv2);

	vector<vector<Ctxt>> convOut2;
	convOperationParallel(stridex, stridey, filterSize, width, length,
					poolOut1[0], noOfFilters, filters, zero, ref(convOut2));
	for (unsigned int i = 1; i < poolOut1.size(); i++) {
		convOperationParallelAdd(stridex, stridey, filterSize, width, length,
				poolOut1[i], noOfFilters, filters, zero, ref(convOut2));
	}
	poolOut1.clear();
	FHE_NTIMER_STOP(conv2);
	tmp = getTimerByName("conv2");
	times = tmp->getTime();
	cout << "Time for conv = " << times << endl;
	cout << "Size after conv: " << convOut2[0].size() << endl;
	cout << "Depth after conv: " << convOut2.size() << endl;

//	for (unsigned int i=0; i < convOut2[1].size(); i++){
//		vector<Ctxt> filteredImage = convOut2[1];
//		for (unsigned int i = 0; i < 5; i++) {
//			ZZX pt;
//			secretKey.Decrypt(pt, filteredImage[i]);
//			cout << pt[0] << " " << flush;
//		}
//		cout << endl;
//	}

	cout << "Pool start ... " << flush;

	width = 11;
	length = 11;

	FHE_NTIMER_START(pool2);

	vector<vector<Ctxt>> poolOut2;
	//poolOperationParallel(poolSize, width, length, zero, convOut1, ref(poolOut1));

	for (unsigned int i = 0; i < convOut2.size(); i++) {
		vector<Ctxt> pooledImage;
		poolOut2.push_back(pooledImage);
	}
	tt = new thread[convOut2.size()];
	for (unsigned int i = 0; i < convOut2.size(); i++) {
		tt[i] = thread(poolEncryptedImage, poolSize, width, length, convOut2[i], zero, ref(poolOut2), i);
	}
	for(unsigned int i=0; i <convOut2.size(); i++){
		tt[i].join();
	}
	convOut2.clear();
	FHE_NTIMER_STOP(pool2);
	tmp = getTimerByName("pool2");
	times = tmp->getTime();
	cout << "end" << endl;
	cout << "# Time for Pool: " << times << endl;
	cout << "# Depth after pool: " << poolOut2.size() << endl;
	cout << "# Size after pool: " << poolOut2[0].size() << endl;
	return 0;
}
