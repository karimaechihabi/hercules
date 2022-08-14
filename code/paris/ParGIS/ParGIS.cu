#include <float.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdbool.h>
#include <cuda.h>
#include "cuda_runtime.h"
#include "ParGIS.h"
__constant__ float sax_breakpointingpu[] ={-2.10003938216135,-2.06623181534436,-2.03283385699339,-1.99984225624759,-1.96725376224613,-1.93506512412816,-1.90327309103284,-1.87187441209932,-1.84086583646678,-1.81024411327436,-1.78000599166123,-1.75014822076655,-1.72066754972946,-1.69156072768914,-1.66282450378474,-1.63445562715541,-1.60645084694032,-1.57880691227863,-1.55152057230949,-1.52458857617206,-1.49800767300551,-1.47177461194898,-1.44588614214165,-1.42033901272266,-1.39512997283118,-1.37025577160636,-1.34571315818736,-1.32149888171335,-1.29760969132347,-1.27404233615690,-1.25079356535278,-1.22786012805028,-1.20523877338855,-1.18292625050675,-1.16091930854405,-1.13921469663960,-1.11780916393255,-1.09669945956208,-1.07588233266733,-1.05535453238746,-1.03511280786164,-1.01515390822902,-0.995474582628758,-0.976071580200017,-0.956941650081954,-0.938081541413728,-0.919488003334496,-0.901157784983417,-0.883087635499651,-0.865274304022356,-0.847714539690689,-0.830405091643811,-0.813342709020879,-0.796524140961053,-0.779946136603490,-0.763605445087349,-0.747498815551789,-0.731622997135969,-0.715974738979046,-0.700550790220181,-0.685347899998530,-0.670362817453254,-0.655592291723510,-0.641033071948457,-0.626681907267254,-0.612535546819059,-0.598590739743031,-0.584844235178328,-0.571292782264110,-0.557933130139534,-0.544762027943759,-0.531776224815945,-0.518972469895249,-0.506347512320829,-0.493898101231846,-0.481620985767457,-0.469512915066820,-0.457570638269096,-0.445790904513441,-0.434170462939015,-0.422706062684976,-0.411394452890483,-0.400232382694694,-0.389216601236769,-0.378343857655865,-0.367610901091141,-0.357014480681756,-0.346551345566869,-0.336218244885638,-0.326011927777221,-0.315929143380778,-0.305966640835466,-0.296121169280445,-0.286389477854873,-0.276768315697909,-0.267254431948711,-0.257844575746437,-0.248535496230247,-0.239323942539300,-0.230206663812752,-0.221180409189764,-0.212241927809494,-0.203387968811101,-0.194615281333742,-0.185920614516577,-0.177300717498764,-0.168752339419462,-0.160272229417829,-0.151857136633025,-0.143503810204207,-0.135208999270534,-0.126969452971165,-0.118781920445259,-0.110643150831974,-0.102549893270468,-0.0944988968999002,-0.0864869108594294,-0.0785106842882140,-0.0705669663254127,-0.0626525061101840,-0.0547640527816864,-0.0468983554790787,-0.0390521633415194,-0.0312222255081671,-0.0234052911181804,-0.0155981093107179,-0.00779742922493818,0.00,0.00779742922493841,0.0155981093107181,0.0234052911181806,0.0312222255081673,0.0390521633415196,0.0468983554790789,0.0547640527816867,0.0626525061101842,0.0705669663254129,0.0785106842882143,0.0864869108594296,0.0944988968999004,0.102549893270468,0.110643150831974,0.118781920445259,0.126969452971166,0.135208999270535,0.143503810204207,0.151857136633025,0.160272229417830,0.168752339419462,0.177300717498764,0.185920614516577,0.194615281333742,0.203387968811101,0.212241927809494,0.221180409189765,0.230206663812752,0.239323942539300,0.248535496230247,0.257844575746437,0.267254431948711,0.276768315697909,0.286389477854873,0.296121169280445,0.305966640835466,0.315929143380778,0.326011927777221,0.336218244885638,0.346551345566869,0.357014480681756,0.367610901091141,0.378343857655865,0.389216601236769,0.400232382694694,0.411394452890483,0.422706062684976,0.434170462939015,0.445790904513441,0.457570638269096,0.469512915066821,0.481620985767457,0.493898101231846,0.506347512320829,0.518972469895249,0.531776224815945,0.544762027943760,0.557933130139534,0.571292782264110,0.584844235178328,0.598590739743031,0.612535546819059,0.626681907267254,0.641033071948457,0.655592291723510,0.670362817453254,0.685347899998530,0.700550790220181,0.715974738979046,0.731622997135969,0.747498815551789,0.763605445087349,0.779946136603489,0.796524140961053,0.813342709020879,0.830405091643811,0.847714539690689,0.865274304022356,0.883087635499651,0.901157784983417,0.919488003334496,0.938081541413728,0.956941650081954,0.976071580200017,0.995474582628758,1.01515390822902,1.03511280786164,1.05535453238746,1.07588233266733,1.09669945956208,1.11780916393255,1.13921469663960,1.16091930854405,1.18292625050675,1.20523877338855,1.22786012805028,1.25079356535278,1.27404233615690,1.29760969132347,1.32149888171335,1.34571315818736,1.37025577160636,1.39512997283118,1.42033901272266,1.44588614214165,1.47177461194898,1.49800767300551,1.52458857617206,1.55152057230949,1.57880691227863,1.60645084694032,1.63445562715541,1.66282450378474,1.69156072768914,1.72066754972946,1.75014822076655,1.78000599166123,1.81024411327436,1.84086583646678,1.87187441209932,1.90327309103284,1.93506512412816,1.96725376224613,1.99984225624759,2.03283385699339,2.06623181534436,2.10003938216135};

__global__ void calculate_lbdold(const sax_type * const saxarray,const float * const paa, const long int M, const int N,float * const sax_breakpoints,bool * positionarray,const float BSF,float segmentsize) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;


	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -


	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) 
		{
			sax_type v = saxarray[j*N+i];
			sax_type region_lower = v;//shift operation 
			sax_type region_upper = (~((int)MAXFLOAT) | region_lower);
			if (region_lower == 0)
			{
				breakpoint_lower = -2000000;
			}
			else
			{
				breakpoint_lower = sax_breakpointingpu[region_lower];//(float)(region_lower-128)*(region_lower-128)/16484.0f;//sax_breakpoints[region_lower];
			}

			if (region_upper == 256 - 1) 
			{
				breakpoint_upper = +2000000;
			}
			else
			{
				breakpoint_upper = sax_breakpointingpu[region_lower+1];//(float)(region_upper+1-128)*(region_upper+1-128)/16484.0f;//sax_breakpoints[region_upper+1];//search in a list(why?)
			}


			if (breakpoint_lower > paa[i]) 
			{
				distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
			}
			else if(breakpoint_upper < paa[i])
			{
				distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
			}

		}

		if(segmentsize*distance<BSF)
		{positionarray[j]=true;}
		else
		{positionarray[j]=false;}
	}
}

__global__ void calculate_lbd(const sax_type * const saxarray,const float * const paa, const long int M, const int N,bool * positionarray,const float BSF,float segmentsize) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;
	float lbsf=BSF/segmentsize;

	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -

	
	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) {
                	if(distance<lbsf)
		{
        	
        		sax_type v = saxarray[j*N+i];

        		sax_type region_lower = v ;//shift operation 
        		sax_type region_upper = (~((int)MAXFLOAT) | region_lower);



        	
        		if (region_lower == 0)
			{
	            		breakpoint_lower = -2000000;
				float breaku=((float)region_lower-127.0f)/128.0f;
            			breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
				if(breakpoint_upper < paa[i])
				{
            				distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
        			}
			}
        		else if (region_upper == 256 - 1) 
			{
            			breakpoint_upper = +2000000;
				float breakx=((float)region_lower-128.0f)/128.0f;
           			breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
				if (breakpoint_lower > paa[i]) 
				{
            				distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
        			}
        		}
        		else 
			{
				float breakx=((float)region_lower-128.0f)/128.0f;
           			breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
				if (breakpoint_lower > paa[i]) 
				{
            				distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
        			}
				else
				{
					float breaku=((float)region_lower-127.0f)/128.0f;
            				breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
					if(breakpoint_upper < paa[i])
					{
            					distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
        				}
        			} 

        		}
		}
						        		
    	}

		if(distance<lbsf)
		{positionarray[j]=true;}
		else
		{positionarray[j]=false;}
	}
}



__global__ void calculate_lbdfloat(const sax_type * const saxarray,const float * const paa, const long int M, const int N,float * positionarray,const float BSF, float segmentsize) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;
	

	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -


	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) 
		{
			if(segmentsize*distance<BSF)
			{
				sax_type v = saxarray[j*N+i];
				sax_type region_lower = v ;//shift operation 
				sax_type region_upper = (~((int)MAXFLOAT) | region_lower);
				if (region_lower == 0)
				{
					breakpoint_lower = -2000000;
					float breaku=((float)region_lower-127.0f)/128.0f;
					breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
					if(breakpoint_upper < paa[i])
					{
						distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
					}
				}
				else if (region_upper == 256 - 1)
				{
					breakpoint_upper = +2000000;
					float breakx=((float)region_lower-128.0f)/128.0f;
					breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
					if (breakpoint_lower > paa[i]) 
					{
						distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
					}
				}
				else
				{
					float breakx=((float)region_lower-128.0f)/128.0f;
					breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
					if (breakpoint_lower > paa[i]) 
					{
						distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
					}
					else
					{
						float breaku=((float)region_lower-127.0f)/128.0f;
						breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
						if(breakpoint_upper < paa[i])
						{
							distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
						}
					}
				}
			}
		}
		positionarray[j]=segmentsize*distance;
	}
}


__global__ void calculate_lbdfloattable(const sax_type * const saxarray,const float * const paa, const long int M, const int N,float * positionarray,const float BSF, float segmentsize) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;
	

	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -


	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) 
		{
			if(segmentsize*distance<BSF)
			{
				sax_type v = saxarray[j*N+i];
				sax_type region_lower = v ;//shift operation 
				sax_type region_upper = (~((int)MAXFLOAT) | region_lower);
				if (region_lower == 0)
				{
					breakpoint_lower = -2000000;
					breakpoint_upper = sax_breakpointingpu[region_upper+1];//search in a list(why?)
					if(breakpoint_upper < paa[i])
					{
						distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
					}
				}
				else if (region_upper == 256 - 1)
				{
					breakpoint_upper = +2000000;
					breakpoint_lower = sax_breakpointingpu[region_lower];
					if (breakpoint_lower > paa[i]) 
					{
						distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
					}
				}
				else
				{
					breakpoint_lower = sax_breakpointingpu[region_lower];
					if (breakpoint_lower > paa[i]) 
					{
						distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
					}
					else
					{
						breakpoint_upper = sax_breakpointingpu[region_upper+1];//search in a list(why?)
						if(breakpoint_upper < paa[i])
						{
							distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
						}
					}
				}
			}
		}
		positionarray[j]=segmentsize*distance;
	}
}






extern "C" float* initialGPU(float *qts, float *gqts, sax_type *saxarray, sax_type *gsaxarray, float *dictionary, float *gdictionary,unsigned long datasize,float *sax_breakpoints )
{

	cudaSetDevice(0); 
	cudaMalloc(&gqts, sizeof(float)*16); 
	return gqts;
}
extern "C" float* initialgqts(float *gqts)
{
		cudaMalloc(&gqts, sizeof(float)*16); 
		return gqts;
}
extern "C" void GPUsyn()
{
cudaDeviceSynchronize();
}




extern "C" float* initialgdictionary(float *gdictionary)
{
		cudaMalloc(&gdictionary, sizeof(float)*FULLSIZE);

		return gdictionary;
}
extern "C" bool* initialgposbitmap(bool *gposbitmap,unsigned long datasize)
{
	cudaMalloc(&gposbitmap, sizeof(bool)*datasize); 
		return gposbitmap;
}
extern "C" bool* initialposbitmap(bool *posbitmap,unsigned long datasize)
{
	cudaMallocHost(&posbitmap, sizeof(bool)*datasize); 
		return posbitmap;
}
extern "C" float* initialgposbitmapfloat(float *gposbitmap,unsigned long datasize)
{
	cudaMalloc(&gposbitmap, sizeof(float)*datasize); 
		return gposbitmap;
}
extern "C" float* initialposbitmapfloat(float *posbitmap,unsigned long datasize)
{
cudaMallocHost(&posbitmap, sizeof(float)*datasize); 
		return posbitmap;
}

extern "C" sax_type* initialgsaxarray(sax_type *gsaxarray,unsigned long datasize)
{
	cudaMalloc(&gsaxarray, sizeof(sax_type)*datasize*16); 

		return gsaxarray;
}
extern "C" sax_type* initialsaxarray(sax_type *saxarray,unsigned long datasize)
{
	cudaMallocHost(&saxarray, sizeof(sax_type)*datasize*16); 

		return saxarray;
}
extern "C" void initialdevice()
{
	cudaSetDevice(0);
}
extern "C" void gpumemcpy(sax_type *gsaxarray,sax_type *saxarray,unsigned long datasize)
{	

	cudaMemcpy(gsaxarray, saxarray,sizeof(sax_type)*datasize*16,cudaMemcpyHostToDevice);
}
extern "C" void gpusaxgridmemcpy(sax_type *gsaxarray,sax_type *saxarray,unsigned long datasize)
{	
	cudaMemcpy(gsaxarray, saxarray,sizeof(sax_type)*datasize*16,cudaMemcpyHostToDevice);
}

extern "C" void gpudictionarymemcpy(float *gdictionary,float *sax_breakpoints)
{	
int offset = ((256 - 1) * (256 - 2)) / 2;
	cudaMemcpy(gdictionary, &sax_breakpoints[offset-1], sizeof(float)*FULLSIZE,cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}





extern "C" void GPUfree(void *devicememorypointer)
{
	cudaFree(devicememorypointer);
}






extern "C" void LBDfloatstreamGPU(sax_type *saxarray, float *posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,float * gposbitmap,int segmentnumber,float segmentsize)
{
	int streamnumber=10;
	cudaMemcpy(gqts, qts,sizeof(float)*segmentnumber,cudaMemcpyHostToDevice);
	cudaStream_t streams[streamnumber];
	for(int i=0;i<streamnumber;i++)
	{
		cudaStreamCreate(&streams[i]);
	}

	for(int i=0;i<streamnumber;i++)
	{
		calculate_lbdfloattable<<<200,500,10,streams[i]>>> (saxarray+i*datasize*segmentnumber/streamnumber,gqts, datasize/streamnumber, segmentnumber, gposbitmap+i*datasize/streamnumber,BSF,segmentsize); 
        cudaMemcpyAsync(posbitmap+i*datasize/streamnumber, gposbitmap+i*datasize/streamnumber, sizeof(float)*datasize/streamnumber,cudaMemcpyDeviceToHost,streams[i]);
	}
	cudaDeviceSynchronize();
}




extern "C" void LBDstreamGPU(sax_type *saxarray, bool *posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,bool * gposbitmap,int segmentnumber,float segmentsize)
{
	int streamnumber=20;
	cudaMemcpy(gqts, qts,sizeof(float)*segmentnumber,cudaMemcpyHostToDevice);
	cudaStream_t streams[streamnumber];
	for(int i=0;i<streamnumber;i++)
	{
		cudaStreamCreate(&streams[i]);
	}
	for(int i=0;i<streamnumber;i++)
	{
		calculate_lbd<<<200,500,10,streams[i]>>> (saxarray+i*datasize*segmentnumber/streamnumber,gqts, datasize/streamnumber, segmentnumber,gposbitmap+i*datasize/streamnumber,BSF, segmentsize); 
		cudaMemcpyAsync(posbitmap+i*datasize/streamnumber, gposbitmap+i*datasize/streamnumber, sizeof(bool)*datasize/streamnumber,cudaMemcpyDeviceToHost,streams[i]);
	}
	cudaDeviceSynchronize();
}


extern "C" void SIMSlowertableGPU(sax_type *saxarray, bool *posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,bool * gposbitmap,int segmentnumber,float segmentsize,float *gdictionary)
{
	int streamnumber=20;
	cudaMemcpy(gqts, qts,sizeof(float)*segmentnumber,cudaMemcpyHostToDevice);
	cudaStream_t streams[streamnumber];
	for(int i=0;i<streamnumber;i++)
	{
		cudaStreamCreate(&streams[i]);
	}
	for(int i=0;i<streamnumber;i++)
	{
		//calculate_lbd<<<200,500,10,streams[i]>>> (saxarray+i*datasize*segmentnumber/streamnumber,gqts, datasize/streamnumber, segmentnumber,gposbitmap+i*datasize/streamnumber,BSF, segmentsize); 
		calculate_lbdold<<<200,500,10,streams[i]>>>(saxarray+i*datasize*segmentnumber/streamnumber,gqts, datasize/streamnumber, segmentnumber,gdictionary,gposbitmap+i*datasize/streamnumber,BSF, segmentsize);
		cudaMemcpyAsync(posbitmap+i*datasize/streamnumber, gposbitmap+i*datasize/streamnumber, sizeof(bool)*datasize/streamnumber,cudaMemcpyDeviceToHost,streams[i]);
	}
	cudaDeviceSynchronize();
}

















/*
extern "C" void SIMSlowerGPU(sax_type *saxarray, bool *posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,bool * gposbitmap,float *gdictionary)
{

	cudaMemcpy(gqts, qts,sizeof(float)*16,cudaMemcpyHostToDevice);
	calculate_lbd<<<200,200,10>>> (saxarray,gqts, datasize, 16, gposbitmap,BSF,16); 

	cudaMemcpy(posbitmap, gposbitmap, sizeof(bool)*datasize,cudaMemcpyDeviceToHost);
//cudaDeviceSynchronize();


}
extern "C" void SIMSlowerGPUfloat(sax_type *saxarray, float *posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,float * gposbitmap,float *gdictionary)
{

	cudaMemcpy(gqts, qts,sizeof(float)*16,cudaMemcpyHostToDevice);
	calculate_lbdfloat<<<200,500,10>>> (saxarray,gqts, datasize, 16,gdictionary, gposbitmap,BSF,16); 
	cudaMemcpy(posbitmap, gposbitmap, sizeof(float)*datasize,cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

extern "C" void initialGPU_grid2(float *qts, float *gqts, sax_type **saxarray, sax_type **gsaxarray, bool **posbitmap, bool **gposbitmap, float *dictionary, float *gdictionary, unsigned long datasize,unsigned long *gridsize,float *sax_breakpoints )
{
	cudaSetDevice(0); 

	for(int i =0;i<65536;i++)
	{
		cudaMallocHost(&saxarray[i], sizeof(sax_type)*gridsize[i]*16); 
		cudaMalloc(&gposbitmap[i], sizeof(float)*gridsize[i]); 
		cudaMalloc(&gsaxarray[i], sizeof(float)*gridsize[i]*16); 
		cudaMallocHost(&posbitmap[i], sizeof(float)*gridsize[i]); 
		cudaMemcpy(gsaxarray[i], saxarray[i],sizeof(float)*gridsize[i]*16,cudaMemcpyHostToDevice);
	}
	cudaMallocHost(&dictionary, 257*sizeof(float)); 
	cudaMallocHost(&qts, sizeof(float)*256); 
	int offset = ((256 - 1) * (256 - 2)) / 2;
	//memcpy(dictionary,&sax_breakpoints[offset-1],sizeof(float)*257);

	cudaMalloc(&gdictionary, sizeof(float)*257);
	cudaMalloc(&gqts, sizeof(float)*256); 

	cudaMemcpy(gdictionary, &sax_breakpoints[offset-1], sizeof(float)*257,cudaMemcpyHostToDevice);
}
extern "C" void SIMSlowerGPUgridstream(sax_type **saxarray, bool **posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,bool ** gposbitmap,float *gdictionary,long int *nodemap, unsigned long *gridnumber)
{
	int streamnumber=(int)datasize;
	cudaStream_t streams[streamnumber];
	cudaMemcpy(gqts, qts,sizeof(float)*16,cudaMemcpyHostToDevice);
	for(int i=0;i<streamnumber;i++)
	{
		cudaStreamCreate(&streams[i]);
	}

	for(long int i=0;i<streamnumber;i++)
	{

		calculate_lbd2<<<500,100,10,streams[i]>>> (saxarray[nodemap[i]],gqts, gridnumber[nodemap[i]], 16,gdictionary, gposbitmap[nodemap[i]],BSF); 
        cudaMemcpyAsync(posbitmap[nodemap[i]], gposbitmap[nodemap[i]], sizeof(bool)*gridnumber[nodemap[i]],cudaMemcpyDeviceToHost,streams[i]);
	}
	for(long int i=0;i<datasize;i++)
	{
		//cudaMemcpy(posbitmap[nodemap[i]], gposbitmap[nodemap[i]], sizeof(bool)*gridnumber[nodemap[i]],cudaMemcpyDeviceToHost,streams[i]);
	}
	cudaDeviceSynchronize();
}
extern "C" void SIMSstreamlowerGPU2(sax_type *saxarray, bool *posbitmap,float * qts,float * gqts,float BSF,unsigned long int startnumber,unsigned long int stopnomber,bool * gposbitmap,float *gdictionary)
{

	int streamnumber=10;
	cudaMemcpy(gqts, qts,sizeof(float)*16,cudaMemcpyHostToDevice);
	unsigned long int datasize=100000000;

	int aaa=startnumber/100000000;
	int bbb=stopnomber/100000000+1;

	if(bbb>=10)
	bbb=10;


	cudaStream_t streams[streamnumber];
	for(int i=0;i<streamnumber;i++)
	{
		cudaStreamCreate(&streams[i]);
	}
	for(int i=0;i<streamnumber;i++)
	{
		calculate_lbd3<<<200,500,10,streams[i]>>> (saxarray+i*datasize*16/streamnumber,gqts, datasize/streamnumber, 16,gposbitmap+i*datasize/streamnumber,BSF,16.0); 
       	cudaMemcpyAsync(posbitmap+i*datasize/streamnumber, gposbitmap+i*datasize/streamnumber, sizeof(bool)*datasize/streamnumber,cudaMemcpyDeviceToHost,streams[i]);
	}

	cudaDeviceSynchronize();

}
extern "C" void SIMSlowerGPUgrid(sax_type **saxarray, bool **posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,bool ** gposbitmap,float *gdictionary,long int *nodemap, unsigned long *gridnumber)
{
	cudaMemcpy(gqts, qts,sizeof(float)*16,cudaMemcpyHostToDevice);
	for(long int i=0;i<datasize;i++)
	{

		calculate_lbd2<<<500,200,10>>> (saxarray[nodemap[i]],gqts, gridnumber[nodemap[i]], 16,gdictionary, gposbitmap[nodemap[i]],BSF); 
	}
	for(long int i=0;i<datasize;i++)
	{
		cudaMemcpy(posbitmap[nodemap[i]], gposbitmap[nodemap[i]], sizeof(bool)*gridnumber[nodemap[i]],cudaMemcpyDeviceToHost);
	}
	cudaDeviceSynchronize();
}

__global__ void calculate_lbd4(const sax_type * const saxarray,const float * const paa, const long int M, const int N,bool * positionarray,const float BSF,const long int offset) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;
	float lbsf=BSF/16.0;

	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -

	
	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) {
                	if(distance<lbsf)
		{
        	
        		sax_type v = saxarray[j*N+i];

        		sax_type region_lower = v ;//shift operation 
        		sax_type region_upper = (~((int)MAXFLOAT) | region_lower);



        	
        		if (region_lower == 0)
			{
	            		breakpoint_lower = -2000000;
				float breaku=((float)region_lower-127.0f)/128.0f;
            			breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
				if(breakpoint_upper < paa[i])
				{
            				distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
        			}
			}
        		else if (region_upper == 256 - 1) 
			{
            			breakpoint_upper = +2000000;
				float breakx=((float)region_lower-128.0f)/128.0f;
           			breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
				if (breakpoint_lower > paa[i]) 
				{
            				distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
        			}
        		}
        		else 
			{
				float breakx=((float)region_lower-128.0f)/128.0f;
           			breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
				if (breakpoint_lower > paa[i]) 
				{
            				distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
        			}
				else
				{
					float breaku=((float)region_lower-127.0f)/128.0f;
            				breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
					if(breakpoint_upper < paa[i])
					{
            					distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
        				}
        			} 

        		}
		}
						        		
    	}

		if(distance<lbsf)
		{positionarray[j]=true;}
		else
		{positionarray[j]=false;}
	}
}


__global__ void calculate_lbd2float(const sax_type * const saxarray,const float * const paa, const long int M, const int N,float * const sax_breakpoints,float * positionarray,const float BSF) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;
	

	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -


	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) 
		{
        	if(16*distance<BSF)
		{
        	
        		sax_type v = saxarray[j*N+i];

        		sax_type region_lower = v ;//shift operation 
        		sax_type region_upper = (~((int)MAXFLOAT) | region_lower);



        	
        		if (region_lower == 0)
			{
	            		breakpoint_lower = -2000000;
        		}
        		else
        		{
				float breakx=((float)region_lower-128.0f)/128.0f;
           			breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
        		}




        		if (region_upper == 256 - 1) 
			{
            			breakpoint_upper = +2000000;
        		}
        		else 
			{
			float breaku=((float)region_lower-127.0f)/128.0f;
            			breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
        		} 

	        			

        		if (breakpoint_lower > paa[i]) 
			{

            			distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
        		}
        		else if(breakpoint_upper < paa[i])
			{
            			distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
        		}
		}
			        				        		
    		}


		positionarray[j]=16*distance;
	}
}

__global__ void calculate_ed2(const float * const Subject, const float * const Query, const long int M, const int N,float * gposbitmap) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float dist=0;
	long int poss;
	//printf("hello jfeowjfw %d\n",thid);
	for (int i = thid; i < M; i += gridDim.x*blockDim.x)
	{	
		dist=0;
		poss=i*N;
		for(int j =0;j<N;j++)
		{
        		dist += (Subject[poss+j]-Query[j])*(Subject[poss+j]-Query[j]);
				
		}
		gposbitmap[i]=dist;
		//if(dist<BSF)
		//resultmap[i]=dist;
	}		
}
__global__ void calculate_lbdold(const sax_type * const saxarray,const float * const paa, const long int M, const int N,float * const sax_breakpoints,bool * positionarray,const float BSF) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;
	

	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -


	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) {
        
        	
        		sax_type v = saxarray[j*N+i];

        		sax_type region_lower = v ;//shift operation 
        		sax_type region_upper = (~((int)MAXFLOAT) | region_lower);



        	
        		if (region_lower == 0)
			{
	            		breakpoint_lower = -2000000;
        		}
        		else
        		{
           			breakpoint_lower = sax_breakpoints[region_lower];//(float)(region_lower-128)*(region_lower-128)/16484.0f;//sax_breakpoints[region_lower];
        		}

        		if (region_upper == 256 - 1) 
			{
            			breakpoint_upper = +2000000;
        		}
        		else
        		{
            			breakpoint_upper = sax_breakpoints[region_lower+1];//(float)(region_upper+1-128)*(region_upper+1-128)/16484.0f;//sax_breakpoints[region_upper+1];//search in a list(why?)
        		}

	        			

        		if (breakpoint_lower > paa[i]) 
			{

            			distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
        		}
        		else if(breakpoint_upper < paa[i])
			{
            			distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
        		}

    		}

		if(16*distance<BSF)
		{positionarray[j]=true;}
		else
		{positionarray[j]=false;}
	}
}
__global__ void calculate_lbd2(const sax_type * const saxarray,const float * const paa, const long int M, const int N,float * const sax_breakpoints,bool * positionarray,const float BSF) 
{
	const int thid = blockDim.x*blockIdx.x + threadIdx.x;
	float distance = 0;
	

	int i=0;
        		float breakpoint_lower = 0; // <-- TODO: calculate breakpoints.
        		float breakpoint_upper = 0; // <-- - || -


	for (int j = thid; j < M; j += gridDim.x*blockDim.x)
	{
		distance=0;
		for (i=0; i<N; i++) {
                	if(16*distance<BSF)
		{
        	
        		sax_type v = saxarray[j*N+i];

        		sax_type region_lower = v ;//shift operation 
        		sax_type region_upper = (~((int)MAXFLOAT) | region_lower);



        	
        		if (region_lower == 0)
			{
	            		breakpoint_lower = -2000000;
        		}
        		else
        		{
				float breakx=((float)region_lower-128.0f)/128.0f;
           			breakpoint_lower = breakx*(breakx*breakx*1.1362582192+0.99800);//sax_breakpoints[region_lower];
        		}

        		if (region_upper == 256 - 1) 
			{
            			breakpoint_upper = +2000000;
        		}
        		else 
			{
			float breaku=((float)region_lower-127.0f)/128.0f;
            			breakpoint_upper = breaku*(1.1362582192*breaku*breaku+0.99800);//sax_breakpoints[region_upper+1];//search in a list(why?)
        		} 

	        			

        		if (breakpoint_lower > paa[i]) 
			{

            			distance += (breakpoint_lower - paa[i])*(breakpoint_lower - paa[i]);
        		}
        		else if(breakpoint_upper < paa[i])
			{
            			distance += (breakpoint_upper - paa[i])*(breakpoint_upper - paa[i]);
        		}
}
						        		
    		}

		if(16*distance<BSF)
		{positionarray[j]=true;}
		else
		{positionarray[j]=false;}
	}
}
extern "C" void SIMSlowerGPUsmall(sax_type *saxarray, bool *posbitmap,float * qts,float * gqts,float BSF,long unsigned datasize,bool * gposbitmap,float *gdictionary)
{

	calculate_lbd<<<20,20,10>>> (saxarray,gqts, datasize, 16,gdictionary, gposbitmap,BSF); 

	cudaMemcpy(posbitmap, gposbitmap, sizeof(bool)*datasize,cudaMemcpyDeviceToHost);

}
extern "C" void copyqts(float * qts,float * gqts)
{
	cudaMemcpy(gqts, qts,sizeof(float)*16,cudaMemcpyHostToDevice);
}
*/