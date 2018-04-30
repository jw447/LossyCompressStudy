/**
 *  @file sz_double_pwr.c
 *  @author Sheng Di and Dingwen Tao
 *  @date Aug, 2016
 *  @brief SZ_Init, Compression and Decompression functions
 * This file contains the compression/decompression functions related to point-wise relative errors
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "sz.h"
#include "CompressElement.h"
#include "DynamicByteArray.h"
#include "DynamicIntArray.h"
#include "TightDataPointStorageD.h"
#include "zlib.h"
#include "rw.h"

void compute_segment_precisions_double_1D(double *oriData, int dataLength,double pwr_err, double* pwrErrBound, unsigned char* pwrErrBoundBytes)
{
	int i = 0, j = 0, k = 0;
	// printf("pw_relBoundRatio=%f\n",pwr_err); // used to be pw_relBoundRatio
	double realPrecision = oriData[0]!=0?fabs(pwr_err*oriData[0]):pwr_err; //initialize a precision = origin[0] * pwr_err
	double approxPrecision;
	unsigned char realPrecBytes[8];
	double curPrecision;
	double curValue;
	for(i=0;i<dataLength;i++) // for all data
	{
		curValue = oriData[i];
		if(i%segment_size==0&&i>0)  // if current value is the first of current segmentation
		{
			//get two first bytes of therealPrecision realPrecision
			doubleToBytes(realPrecBytes, realPrecision);
			// why replace?
			memset(&realPrecBytes[2], 0, 6); // replace with 0 from [2], replace 6 of them to the end
			approxPrecision = bytesToDouble(realPrecBytes);
			//put the realPrecision in float* pwrErBound
			pwrErrBound[j++] = approxPrecision; // this is float number instead of double
			// printf("%.15f\n",pwrErrBound[j]);
			//put the two bytes in pwrErrBoundBytes
			pwrErrBoundBytes[k++] = realPrecBytes[0];
			pwrErrBoundBytes[k++] = realPrecBytes[1];
			realPrecision = curValue!=0?fabs(pwr_err*curValue):pwr_err;
		}
		else if(curValue!=0)
		{
			curPrecision = fabs(pwr_err*curValue);
			if(realPrecision>curPrecision)
				realPrecision = curPrecision;
		}
	}
	doubleToBytes(realPrecBytes, realPrecision);
	memset(&realPrecBytes[2], 0, 6);
	approxPrecision = bytesToDouble(realPrecBytes);
	//put the realPrecision in float* pwrErBound
	pwrErrBound[j++] = approxPrecision;
	//put the two bytes in pwrErrBoundBytes
	pwrErrBoundBytes[k++] = realPrecBytes[0];
	pwrErrBoundBytes[k++] = realPrecBytes[1];
	// printf("%.15f\n",pwrErrBound[j]);
	// printf("%d\n",j);
}

// not understanded
unsigned int optimize_intervals_double_1D_pwr(double *oriData, int dataLength,
	double* pwrErrBound)
{
	int i = 0, j = 0;
	double realPrecision = pwrErrBound[j++];
	unsigned long radiusIndex;
	double pred_value = 0, pred_err;
	//printf("maxRangeRadius=%d\n",maxRangeRadius);//32768 decided in sz.c
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int)); // all 0
	// printf("sampleDistance=%d\n",sampleDistance); // 10?
	int totalSampleSize = dataLength/sampleDistance;
	// printf("totalSampleSize=%d\n",totalSampleSize);
	for(i=2;i<dataLength;i++)
	{
		if(i%segment_size==0)
			realPrecision = pwrErrBound[j++];
		if(i%sampleDistance==0)
		{
			pred_value = 2*oriData[i-1] - oriData[i-2];
			//pred_value = oriData[i-1];
			pred_err = fabs(pred_value - oriData[i]);
			radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
			if(radiusIndex>=maxRangeRadius)
				radiusIndex = maxRangeRadius - 1;
			intervals[radiusIndex]++;
		}
	}
	//compute the appropriate number

	int targetCount = (int)(totalSampleSize*predThreshold); // threshold = 0.99

	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	// printf("i=%d\n",i);
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
	return powerOf2;
}

void compute_segment_precisions_double_2D(double *oriData, double* pwrErrBound,
int r1, int r2, int R2, int edgeSize, unsigned char* pwrErrBoundBytes, double Min, double Max)
{
	int i = 0, j = 0, k = 0, p = 0, index = 0, J; //I=-1,J=-1 if they are needed
	double realPrecision;
	double approxPrecision;
	unsigned char realPrecBytes[8];
	double curValue, curAbsValue;
	double* minAbsValues = (double*)malloc(R2*sizeof(double));

	double max = fabs(Min)<fabs(Max)?fabs(Max):fabs(Min); //get the max abs value.
	for(i=0;i<R2;i++)
		minAbsValues[i] = max;
	for(i=0;i<r1;i++)
	{
		for(j=0;j<r2;j++)
		{
			index = i*r2+j;
			curValue = oriData[index];
			if(((i%edgeSize==edgeSize-1 || i==r1-1) &&j%edgeSize==0&&j>0) || (i%edgeSize==0&&j==0&&i>0))
			{
				realPrecision = pw_relBoundRatio*minAbsValues[J];
				doubleToBytes(realPrecBytes, realPrecision);
				memset(&realPrecBytes[2], 0, 6);
				approxPrecision = bytesToDouble(realPrecBytes);
				//put the realPrecision in double* pwrErBound
				pwrErrBound[p++] = approxPrecision;
				//put the two bytes in pwrErrBoundBytes
				pwrErrBoundBytes[k++] = realPrecBytes[0];
				pwrErrBoundBytes[k++] = realPrecBytes[1];
				minAbsValues[J] = max;
			}
			if(j==0)
				J = 0;
			else if(j%edgeSize==0)
				J++;
			if(curValue!=0)
			{
				curAbsValue = fabs(curValue);
				if(minAbsValues[J]>curAbsValue)
					minAbsValues[J] = curAbsValue;
			}
		}
	}
	realPrecision = pw_relBoundRatio*minAbsValues[J];
	doubleToBytes(realPrecBytes, realPrecision);
	realPrecBytes[2] = realPrecBytes[3] = 0;
	approxPrecision = bytesToDouble(realPrecBytes);
	//put the realPrecision in double* pwrErBound
	pwrErrBound[p++] = approxPrecision;
	//put the two bytes in pwrErrBoundBytes
	pwrErrBoundBytes[k++] = realPrecBytes[0];
	pwrErrBoundBytes[k++] = realPrecBytes[1];

	free(minAbsValues);
}

unsigned int optimize_intervals_double_2D_pwr(double *oriData, int r1, int r2, int R2, int edgeSize, double* pwrErrBound)
{
	int i = 0,j = 0, index, I=0, J=0;
	double realPrecision = pwrErrBound[0];
	unsigned long radiusIndex;
	double pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = r1*r2/sampleDistance;
	int ir2;
	for(i=1;i<r1;i++)
	{
		ir2 = i*r2;
		if(i%edgeSize==0)
		{
			I++;
			J = 0;
		}
		for(j=1;j<r2;j++)
		{
			index = ir2+j;
			if(j%edgeSize==0)
				J++;

			if((i+j)%sampleDistance==0)
			{
				realPrecision = pwrErrBound[I*R2+J];
				pred_value = oriData[index-1] + oriData[index-r2] - oriData[index-r2-1];
				pred_err = fabs(pred_value - oriData[index]);
				radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
				if(radiusIndex>=maxRangeRadius)
					radiusIndex = maxRangeRadius - 1;
				intervals[radiusIndex]++;
			}
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("maxRangeRadius = %d, accIntervals=%d, powerOf2=%d\n", maxRangeRadius, accIntervals, powerOf2);
	return powerOf2;
}

void compute_segment_precisions_double_3D(double *oriData, double* pwrErrBound,
int r1, int r2, int r3, int R2, int R3, int edgeSize, unsigned char* pwrErrBoundBytes, double Min, double Max)
{
	int i = 0, j = 0, k = 0, p = 0, q = 0, index = 0, J = 0, K = 0; //I=-1,J=-1 if they are needed
	int r23 = r2*r3, ir, jr;
	double realPrecision;
	double approxPrecision;
	unsigned char realPrecBytes[8];
	double curValue, curAbsValue;

	double** minAbsValues = create2DArray_double(R2, R3);
	double max = fabs(Min)<fabs(Max)?fabs(Max):fabs(Min); //get the max abs value.
	for(i=0;i<R2;i++)
		for(j=0;j<R3;j++)
			minAbsValues[i][j] = max;

	for(i=0;i<r1;i++)
	{
		ir = i*r23;
		if(i%edgeSize==0&&i>0)
		{
			realPrecision = pw_relBoundRatio*minAbsValues[J][K];
			doubleToBytes(realPrecBytes, realPrecision);
			memset(&realPrecBytes[2], 0, 6);
			approxPrecision = bytesToDouble(realPrecBytes);
			//put the realPrecision in float* pwrErBound
			pwrErrBound[p++] = approxPrecision;
			//put the two bytes in pwrErrBoundBytes
			//printf("q=%d, i=%d, j=%d, k=%d\n",q,i,j,k);
			pwrErrBoundBytes[q++] = realPrecBytes[0];
			pwrErrBoundBytes[q++] = realPrecBytes[1];
			minAbsValues[J][K] = max;
		}
		for(j=0;j<r2;j++)
		{
			jr = j*r3;
			if((i%edgeSize==edgeSize-1 || i == r1-1)&&j%edgeSize==0&&j>0)
			{
				realPrecision = pw_relBoundRatio*minAbsValues[J][K];
				doubleToBytes(realPrecBytes, realPrecision);
				memset(&realPrecBytes[2], 0, 6);
				approxPrecision = bytesToDouble(realPrecBytes);
				//put the realPrecision in float* pwrErBound
				pwrErrBound[p++] = approxPrecision;
				//put the two bytes in pwrErrBoundBytes
				//printf("q=%d, i=%d, j=%d, k=%d\n",q,i,j,k);
				pwrErrBoundBytes[q++] = realPrecBytes[0];
				pwrErrBoundBytes[q++] = realPrecBytes[1];
				minAbsValues[J][K] = max;
			}

			if(j==0)
				J = 0;
			else if(j%edgeSize==0)
				J++;

			for(k=0;k<r3;k++)
			{
				index = ir+jr+k;
				curValue = oriData[index];
				if((i%edgeSize==edgeSize-1 || i == r1-1)&&(j%edgeSize==edgeSize-1||j==r2-1)&&k%edgeSize==0&&k>0)
				{
					realPrecision = pw_relBoundRatio*minAbsValues[J][K];
					doubleToBytes(realPrecBytes, realPrecision);
					memset(&realPrecBytes[2], 0, 6);
					approxPrecision = bytesToDouble(realPrecBytes);
					//put the realPrecision in float* pwrErBound
					pwrErrBound[p++] = approxPrecision;
					//put the two bytes in pwrErrBoundBytes
					//printf("q=%d, i=%d, j=%d, k=%d\n",q,i,j,k);
					pwrErrBoundBytes[q++] = realPrecBytes[0];
					pwrErrBoundBytes[q++] = realPrecBytes[1];
					minAbsValues[J][K] = max;
				}

				if(k==0)
					K = 0;
				else if(k%edgeSize==0)
					K++;

				if(curValue!=0)
				{
					curAbsValue = fabs(curValue);
					if(minAbsValues[J][K]>curAbsValue)
						minAbsValues[J][K] = curAbsValue;
				}
			}
		}
	}

	realPrecision = pw_relBoundRatio*minAbsValues[J][K];
	doubleToBytes(realPrecBytes, realPrecision);
	memset(&realPrecBytes[2], 0, 6);
	approxPrecision = bytesToDouble(realPrecBytes);
	//put the realPrecision in float* pwrErBound
	pwrErrBound[p++] = approxPrecision;
	//put the two bytes in pwrErrBoundBytes
	pwrErrBoundBytes[q++] = realPrecBytes[0];
	pwrErrBoundBytes[q++] = realPrecBytes[1];

	free2DArray_double(minAbsValues, R2);
}

unsigned int optimize_intervals_double_3D_pwr(double *oriData, int r1, int r2, int r3, int R2, int R3, int edgeSize, double* pwrErrBound)
{
	int i,j,k, ir,jr,index, I = 0,J=0,K=0;
	double realPrecision = pwrErrBound[0];
	unsigned long radiusIndex;
	int r23=r2*r3;
	int R23 = R2*R3;
	double pred_value = 0, pred_err;
	int *intervals = (int*)malloc(maxRangeRadius*sizeof(int));
	memset(intervals, 0, maxRangeRadius*sizeof(int));
	int totalSampleSize = r1*r2*r3/sampleDistance;
	for(i=1;i<r1;i++)
	{
		ir = i*r23;
		if(i%edgeSize==0)
		{
			I++;
			J = 0;
		}
		for(j=1;j<r2;j++)
		{
			jr = j*r3;
			if(j%edgeSize==0)
			{
				J++;
				K = 0;
			}
			for(k=1;k<r3;k++)
			{
				index = ir+jr+k;
				if(k%edgeSize==0)
					K++;
				if((i+j+k)%sampleDistance==0)
				{
					realPrecision = pwrErrBound[I*R23+J*R2+K];
					pred_value = oriData[index-1] + oriData[index-r3] + oriData[index-r23]
					- oriData[index-1-r23] - oriData[index-r3-1] - oriData[index-r3-r23] + oriData[index-r3-r23-1];
					pred_err = fabs(pred_value - oriData[index]);
					radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
					if(radiusIndex>=maxRangeRadius)
						radiusIndex = maxRangeRadius - 1;
					intervals[radiusIndex]++;
				}
			}
		}
	}
	//compute the appropriate number
	int targetCount = (int)(totalSampleSize*predThreshold);
	int sum = 0;
	for(i=0;i<maxRangeRadius;i++)
	{
		sum += intervals[i];
		if(sum>targetCount)
			break;
	}
	if(i>=maxRangeRadius)
		i = maxRangeRadius-1;
	unsigned int accIntervals = 2*(i+1);
	unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);

	if(powerOf2<32)
		powerOf2 = 32;

	free(intervals);
	//printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
	return powerOf2;
}

void SZ_compress_args_double_NoCkRngeNoGzip_1D_pwr(unsigned char** newByteData, double *oriData, double pwr_err,
int dataLength, int *outSize, double min, double max)
{
	// printf("1,intvR=%d\n",intvRadius);
	SZ_Reset();
	int pwrLength = dataLength%segment_size==0?dataLength/segment_size:dataLength/segment_size+1; // array length of point-wise relative errorbound
	double* pwrErrBound = (double*)malloc(sizeof(double)*pwrLength); // array of point-wise relative errorbound for each segment
	int pwrErrBoundBytes_size = sizeof(unsigned char)*pwrLength*2;
	unsigned char* pwrErrBoundBytes = (unsigned char*)malloc(pwrErrBoundBytes_size);

	// printf("lengthof pwr_error bound:%ld\n",sizeof(pwrErrBound));
	//for(int i = 0;i<8832;i++){
	//	printf("%.10f\n",pwrErrBound[i]);
	//}

	// printf("dataLength=%d\n",dataLength);
	// oridata is the input data, datalength is the number of data points, pwr_err is the required point-wise
	// errorbound, pwrErrBound is the output of errorbound for each segment, pwrErrBoundBytes is the byte[0]/[1] of pwrErrBound
	compute_segment_precisions_double_1D(oriData, dataLength, pwr_err, pwrErrBound, pwrErrBoundBytes);
	//printf("PointCnt=%d\n",dataLength);
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		// printf("optimization=1\n");
		quantization_intervals = optimize_intervals_double_1D_pwr(oriData, dataLength, pwrErrBound);
		updateQuantizationInfo(quantization_intervals);
		//printf("quantization_intervals=%d\n",quantization_intervals);
		// intvCapacity = quant_intervals;
		// intvRadius = quant_intervals/2;
		// printf("quantization_intervals=%d\n",quantization_intervals);
	}
	else{
		quantization_intervals = intvCapacity;
		// printf("quantization_intervals=%d\n",quantization_intervals);
	}

	//clearHuffmanMem();
	// printf("intvCapacity=%d\n",intvCapacity);
	// printf("2,intvR=%d2\n",intvRadius);
	//printf("quantization_intervals=%d\n",quantization_intervals);
	int i = 0, j = 0, reqLength;
	double realPrecision = pwrErrBound[j++];
	double medianValue = 0;
	double radius = fabs(max)<fabs(min)?fabs(min):fabs(max);
	short radExpo = getExponent_double(radius);

	computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);  // required length of? [12, 64]

	int* type = (int*) malloc(dataLength*sizeof(int));
	//type[dataLength]=0;

	double* spaceFillingValue = oriData; //

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	type[0] = 0;
	//printf("%d\n",type[0]);
	unsigned char preDataBytes[4] = {0};
	intToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;
	double last3CmprsData[3] = {0};

	DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
	// printf("sizeof lossy element=%ld\n",sizeof(lce));
	//add the first data
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_double(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);

	//add the second data
	type[1] = 0;
	//printf("%d\n",type[1]);
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleDoubleValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,4);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	listAdd_double(last3CmprsData, vce->data);
	//printf("%.30G\n",last3CmprsData[0]);
	// printf("3,intvR=%d\n",intvRadius);
	int state;
	double lcf, qcf;
	double checkRadius;
	double curData;
	double pred;
	double predAbsErr;
	double min_pred, minErr, minIndex;
	int a = 0;
	checkRadius = (intvCapacity-1)*realPrecision;
	// printf("%f\n",checkRadius);
	double interval = 2*realPrecision;
	int updateReqLength = 0; //a marker: 1 means already updated

	int hit = 0;
	int miss = 2; // first two values are not hit.
	int reqBits = 0;

	// printf("datalength=%d\n",dataLength);
	// printf("qf:%d\n",quantization_intervals);
	// printf("realPrecision:%f\n",realPrecision);
	// printf("checkRadius:%f\n",checkRadius);
	// printf("interval:%f\n",interval);
	// printf("intvR=%d\n",intvRadius);
	// // printf("%d\n",dataLength);
	// printf("-------------------\n");
  // mean of dataset ----------
	// double tmp = 0.0;
	// for(i=0;i<dataLength;i++){
	// 	// tmp += spaceFillingValue[i];
	// 	printf("%f\n",spaceFillingValue[i]);
	// }
	// printf("%2.10f",realPrecision);
	// printf("%f\n",last3CmprsData[0]);
	// printf("%f\n",last3CmprsData[1]);
	// printf("mean of dataset=%f\n",tmp/dataLength);
  // ------------------------------
	int f = 0;
	for(i=2;i<dataLength;i++)
	{

		curData = spaceFillingValue[i];
		// printf("data:%f\n",curData);
		if(i%segment_size==0)
		{
			realPrecision = pwrErrBound[j++];
			// printf("intvCapacity=%d\n",intvCapacity);
			checkRadius = (intvCapacity-1)*realPrecision;
			//checkradius_s += checkRadius;
			// printf("data=%f\n",curData);
			// printf("%f\n",checkRadius);
			interval = 2*realPrecision;
			updateReqLength = 0;
		}

		pred = 2*last3CmprsData[0] - last3CmprsData[1];

		// printf("%.15f\n",pred);
		//pred = last3CmprsData[0];
		// printf("******************************\n");
		predAbsErr = fabs(curData - pred);
		// printf("%.15f\n",curData);
		// printf("%2.15f\n",(pred));
		// printf("%.15f\n",predAbsErr);
		// printf("%.15f\n",checkRadius);
		// printf("pred:%f\n",pred);
		// printf("predAbsErr:%f\n",predAbsErr);
		// printf("checkRadius:%f\n",checkRadius);
		if(predAbsErr<checkRadius)
		{
			// f++;
			// printf("Hit\n");
			// printf("%d\n",type[i]);
			state = (predAbsErr/realPrecision+1)/2;
			// printf("state:%d\n",state);
			// printf("state=%d\n",state);
			// printf("%f\n",(curData-pred));
			hit = hit + 1;
			if(curData>=pred)
			{
				type[i] = intvRadius+state;
				// printf("%d\n",type[i]);
				// printf(">\n");
				pred = pred + state*interval;
				// printf("pred_updated:%f\n",pred);
			}
			else //curData<pred
			{
				type[i] = intvRadius-state;
				// printf("%d\n",type[i]);
				// printf("type[i]:%d\n",type[i]);
				// printf("<\n");
				pred = pred - state*interval;
				// printf("pred_updated:%f\n",pred);

			}
      //printf("%d\n",type[i]);
			// printf("pred_final:%f\n",pred);
			listAdd_double(last3CmprsData, pred);
			// printf("............................\n");
			continue;
		}

		if(updateReqLength==0)
		{
			computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
			reqBytesLength = reqLength/8;
			resiBitsLength = reqLength%8;
			updateReqLength = 1;
		}

		reqBits += reqLength;
		// printf("Miss\n");
		miss = miss + 1;
		type[i] = 0;
		// printf("%d\n",type[i]);
		// printf("+++++++++++++++++++++++++\n");

		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleDoubleValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);

		memcpy(preDataBytes,vce->curBytes,4);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);

		listAdd_double(last3CmprsData, vce->data);
	}//end of for
	// printf("%d\n",hit);
	// printf("%.20f\n",(double)hit/dataLength);
	//printf("average checkRadius = %.15f\n",checkradius_s/(dataLength/32));
	// printf("f=%d\n",f);
	//printf("----------\n");
	// for(i=0;i<dataLength;i++){
  //           if(type[i]!=0){
	// 	printf("%d\n",type[i]);
  //           }
	// }
	//printf("hit: %d\n",hit);
	//printf("miss: %d\n",miss);
	// printf("hitRatio: %f\n",(float)hit/(float)dataLength);
	// printf("------reqBits=%d\n",reqBits);
	// printf("outlierSize?: %d\n",reqBits/8);
	// printf("hit=%d\n", hit);
	// printf("miss=%d\n",miss);
	//printf("outlierSize: %d\n", outlierSize);
	// printf("sizeof lossy element=%ld\n",sizeof(*lce));
	//	char* expSegmentsInBytes;
	//	int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
	int exactDataNum = exactLeadNumArray->size; // number of missing data
	// printf("exactLeadNumArray=%d\n",sizeof(exactLeadNumArray));
	// printf("exactMidByteArray=%d\n",sizeof(exactMidByteArray));
	// printf("resiBitArray=%d\n",sizeof(resiBitArray));

	TightDataPointStorageD* tdps;
	// printf("2\n");
	new_TightDataPointStorageD(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, pwrErrBoundBytes, pwrErrBoundBytes_size, radExpo);
	// printf("-------insight of tdps\n");
	// printf("dataLength=%d\n",dataLength);
	// printf("exactDataNum=%d\n",exactDataNum);
	// //printf("type=%d\n",type);
	// printf("realPrecision=%f\n",realPrecision);
	// printf("medianValue=%f\n",medianValue);
	// printf("reqLength=%d\n",reqLength);
	// printf("quantization_intervals=%d\n",quantization_intervals);
	// //printf("pwrErrBoundBytes=%d\n",pwrErrBoundBytes);
	// printf("pwrErrBoundBytes_size=%d\n",pwrErrBoundBytes_size);
	// printf("radExpo=%d\n",radExpo);
	// printf("-----------------------\n");
	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);

	convertTDPStoFlatBytes_double(tdps, newByteData, outSize);
	// printf("tdps=%d\n",sizeof(tdps));
	// printf("?=%d\n",exactMidByteArray->size);
	// printf("?=%d\n",resiBitArray->size);
	// printf("?=%d\n",resiBitLengthArray->size);
	// printf("?=%d\n",pwrErrBoundBytes_size);

	//printf("/outSize/=%d\n",*outSize);
	int doubleSize=sizeof(double);
	// printf("-----%d\n",doubleSize*dataLength);
	if(*outSize>dataLength*doubleSize)
	{
		// printf("yes\n");
		int k = 0, i;
		tdps->isLossless = 1;
		int totalByteLength = 3 + 4 + 1 + doubleSize*dataLength;
		*newByteData = (unsigned char*)malloc(totalByteLength);

		unsigned char dsLengthBytes[4];
		intToBytes_bigEndian(dsLengthBytes, dataLength);//4
		for (i = 0; i < 3; i++)//3
			(*newByteData)[k++] = versionNumber[i];
		for (i = 0; i < 4; i++)//4
			(*newByteData)[k++] = dsLengthBytes[i];
		(*newByteData)[k++] = 16;	//=00010000

		if(sysEndianType==BIG_ENDIAN_SYSTEM)
			memcpy((*newByteData)+8, oriData, dataLength*doubleSize);
		else
		{
			unsigned char* p = (*newByteData)+8;
			for(i=0;i<dataLength;i++,p+=doubleSize)
				doubleToBytes(p, oriData[i]);
		}
		*outSize = totalByteLength;
	}
	// printf("pwr:%d\n",*outSize);
	free(pwrErrBound);

	free(vce);
	free(lce);
	free_TightDataPointStorageD(tdps);
	free(exactMidByteArray);
}

/**
 *
 * Note: @r1 is high dimension
 * 		 @r2 is low dimension
 * */
void SZ_compress_args_double_NoCkRngeNoGzip_2D_pwr(unsigned char** newByteData, double *oriData, int r1, int r2,
int *outSize, double min, double max)
{
	SZ_Reset();
	int dataLength=r1*r2;
	int blockEdgeSize = computeBlockEdgeSize_2D(segment_size);
	int R1 = 1+(r1-1)/blockEdgeSize;
	int R2 = 1+(r2-1)/blockEdgeSize;
	double* pwrErrBound = (double*)malloc(sizeof(double)*R1*R2);
	int pwrErrBoundBytes_size = sizeof(unsigned char)*R1*R2*2;
	unsigned char* pwrErrBoundBytes = (unsigned char*)malloc(pwrErrBoundBytes_size);

	compute_segment_precisions_double_2D(oriData, pwrErrBound, r1, r2, R2, blockEdgeSize, pwrErrBoundBytes, min, max);
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_double_2D_pwr(oriData, r1, r2, R2, blockEdgeSize, pwrErrBound);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;
	//clearHuffmanMem();
	//printf("quantization_intervals=%d\n",quantization_intervals);

	int i=0,j=0,I=0,J=0,reqLength;
	double realPrecision = pwrErrBound[I*R2+J];
	double pred1D, pred2D;
	double diff = 0.0;
	double itvNum = 0;
	double *P0, *P1;

	P0 = (double*)malloc(r2*sizeof(double));
	memset(P0, 0, r2*sizeof(double));
	P1 = (double*)malloc(r2*sizeof(double));
	memset(P1, 0, r2*sizeof(double));

	double medianValue = 0;
	double radius = fabs(max)<fabs(min)?fabs(min):fabs(max);
	short radExpo = getExponent_double(radius);
	int updateReqLength = 1;

	computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));
	//type[dataLength]=0;

	double* spaceFillingValue = oriData; //

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	type[0] = 0;

	unsigned char preDataBytes[8];
	longToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));

	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,8);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum =  fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
	}
	else
	{
		type[1] = 0;

		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleDoubleValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,8);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r2-1 */
	for (j = 2; j < r2; j++)
	{
		if(j%blockEdgeSize==0)
		{
			J++;
			realPrecision = pwrErrBound[I*R2+J];
			updateReqLength = 0;
		}

		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
		}
		else
		{
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}

			type[j] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r1-1 */
	int index;
	for (i = 1; i < r1; i++)
	{
		/* Process row-i data 0 */
		index = i*r2;
		J = 0;
		if(i%blockEdgeSize==0)
			I++;
		realPrecision = pwrErrBound[I*R2+J]; //J==0
		updateReqLength = 0;

		pred1D = P1[0];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}

			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}

		/* Process row-i data 1 --> r2-1*/
		for (j = 1; j < r2; j++)
		{
			index = i*r2+j;
			if(j%blockEdgeSize==0)
			{
				J++;
				realPrecision = pwrErrBound[I*R2+J];
				updateReqLength = 0;
			}
			pred2D = P0[j-1] + P1[j] - P1[j-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}

				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

		double *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}

	if(r2!=1)
		free(P0);
	free(P1);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageD* tdps;

	new_TightDataPointStorageD(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, pwrErrBoundBytes, pwrErrBoundBytes_size, radExpo);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);

	convertTDPStoFlatBytes_double(tdps, newByteData, outSize);

	free(pwrErrBound);

	free(vce);
	free(lce);
	free_TightDataPointStorageD(tdps);
	free(exactMidByteArray);
}

void SZ_compress_args_double_NoCkRngeNoGzip_3D_pwr(unsigned char** newByteData, double *oriData, int r1, int r2, int r3, int *outSize, double min, double max)
{
	SZ_Reset();
	int dataLength=r1*r2*r3;

	int blockEdgeSize = computeBlockEdgeSize_3D(segment_size);
	int R1 = 1+(r1-1)/blockEdgeSize;
	int R2 = 1+(r2-1)/blockEdgeSize;
	int R3 = 1+(r3-1)/blockEdgeSize;
	double* pwrErrBound = (double*)malloc(sizeof(double)*R1*R2*R3);
	int pwrErrBoundBytes_size = sizeof(unsigned char)*R1*R2*R3*2;
	unsigned char* pwrErrBoundBytes = (unsigned char*)malloc(pwrErrBoundBytes_size);

	compute_segment_precisions_double_3D(oriData, pwrErrBound, r1, r2, r3, R2, R3, blockEdgeSize, pwrErrBoundBytes, min, max);
	unsigned int quantization_intervals;
	if(optQuantMode==1)
	{
		quantization_intervals = optimize_intervals_double_3D_pwr(oriData, r1, r2, r3, R2, R3, blockEdgeSize, pwrErrBound);
		updateQuantizationInfo(quantization_intervals);
	}
	else
		quantization_intervals = intvCapacity;
	//clearHuffmanMem();
	int i=0,j=0,k=0, reqLength, I = 0, J = 0, K = 0;
	double realPrecision = pwrErrBound[0];
	double pred1D, pred2D, pred3D;
	double diff = 0.0;
	double itvNum = 0;
	double *P0, *P1;

	int r23 = r2*r3;
	int R23 = R2*R3;
	P0 = (double*)malloc(r23*sizeof(double));
	P1 = (double*)malloc(r23*sizeof(double));
	double radius = fabs(max)<fabs(min)?fabs(min):fabs(max);
	double medianValue = 0;
	short radExpo = getExponent_double(radius);
	int updateReqLength = 0;
	computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);

	int* type = (int*) malloc(dataLength*sizeof(int));
	//type[dataLength]=0;

	double* spaceFillingValue = oriData; //

	DynamicByteArray *resiBitLengthArray;
	new_DBA(&resiBitLengthArray, DynArrayInitLen);

	DynamicIntArray *exactLeadNumArray;
	new_DIA(&exactLeadNumArray, DynArrayInitLen);

	DynamicByteArray *exactMidByteArray;
	new_DBA(&exactMidByteArray, DynArrayInitLen);

	DynamicIntArray *resiBitArray;
	new_DIA(&resiBitArray, DynArrayInitLen);

	type[0] = 0;

	unsigned char preDataBytes[8];
	longToBytes_bigEndian(preDataBytes, 0);

	int reqBytesLength = reqLength/8;
	int resiBitsLength = reqLength%8;

	DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
	LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));


	///////////////////////////	Process layer-0 ///////////////////////////
	/* Process Row-0 data 0*/
	type[0] = 0;
	addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
	compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
	updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
	memcpy(preDataBytes,vce->curBytes,8);
	addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
	P1[0] = vce->data;

	/* Process Row-0 data 1*/
	pred1D = P1[0];
	diff = spaceFillingValue[1] - pred1D;

	itvNum = fabs(diff)/realPrecision + 1;

	if (itvNum < intvCapacity)
	{
		if (diff < 0) itvNum = -itvNum;
		type[1] = (int) (itvNum/2) + intvRadius;
		P1[1] = pred1D + 2 * (type[1] - intvRadius) * realPrecision;
	}
	else
	{
		if(updateReqLength==0)
		{
			computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
			reqBytesLength = reqLength/8;
			resiBitsLength = reqLength%8;
			updateReqLength = 1;
		}

		type[1] = 0;

		addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
		compressSingleDoubleValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
		updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
		memcpy(preDataBytes,vce->curBytes,8);
		addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
		P1[1] = vce->data;
	}

    /* Process Row-0 data 2 --> data r3-1 */
	for (j = 2; j < r3; j++)
	{
		if(j%blockEdgeSize==0)
		{
			J++;
			realPrecision = pwrErrBound[J];
			updateReqLength = 0;
		}
		pred1D = 2*P1[j-1] - P1[j-2];
		diff = spaceFillingValue[j] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[j] = (int) (itvNum/2) + intvRadius;
			P1[j] = pred1D + 2 * (type[j] - intvRadius) * realPrecision;
		}
		else
		{
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}

			type[j] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[j], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[j] = vce->data;
		}
	}

	/* Process Row-1 --> Row-r2-1 */
	int index;
	K = 0;
	for (i = 1; i < r2; i++)
	{
		/* Process row-i data 0 */
		index = i*r3;

		J = 0;
		if(i%blockEdgeSize==0)
			I++;
		realPrecision = pwrErrBound[I*R3+J]; //J==0
		updateReqLength = 0;

		pred1D = P1[index-r3];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P1[index] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}

			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P1[index] = vce->data;
		}

		/* Process row-i data 1 --> data r3-1*/
		for (j = 1; j < r3; j++) //note that this j refers to fastest dimension (lowest order)
		{
			index = i*r3+j;
			if(j%blockEdgeSize==0)
			{
				J++;
				realPrecision = pwrErrBound[I*R3+J];
				updateReqLength = 0;
			}

			pred2D = P1[index-1] + P1[index-r3] - P1[index-r3-1];

			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P1[index] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}

				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P1[index] = vce->data;
			}
		}
	}


	///////////////////////////	Process layer-1 --> layer-r1-1 ///////////////////////////

	for (k = 1; k < r1; k++)
	{
		/* Process Row-0 data 0*/
		index = k*r23;
		I = 0;
		J = 0;
		if(k%blockEdgeSize==0)
			K++;
		realPrecision = pwrErrBound[K*R23]; //J==0
		updateReqLength = 0;

		pred1D = P1[0];
		diff = spaceFillingValue[index] - pred1D;

		itvNum = fabs(diff)/realPrecision + 1;

		if (itvNum < intvCapacity)
		{
			if (diff < 0) itvNum = -itvNum;
			type[index] = (int) (itvNum/2) + intvRadius;
			P0[0] = pred1D + 2 * (type[index] - intvRadius) * realPrecision;
		}
		else
		{
			if(updateReqLength==0)
			{
				computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
				reqBytesLength = reqLength/8;
				resiBitsLength = reqLength%8;
				updateReqLength = 1;
			}

			type[index] = 0;

			addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
			compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
			updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
			memcpy(preDataBytes,vce->curBytes,8);
			addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
			P0[0] = vce->data;
		}


	    /* Process Row-0 data 1 --> data r3-1 */
		for (j = 1; j < r3; j++)
		{
			index = k*r23+j;
				if(index==46788)
					printf("index=%d\n", index);
			if(j%blockEdgeSize==0)
			{
				J++;
				realPrecision = pwrErrBound[K*R23+J];
				updateReqLength = 0;
			}
			pred2D = P0[j-1] + P1[j] - P1[j-1];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[j] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}

				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[j] = vce->data;
			}
		}

	    /* Process Row-1 --> Row-r2-1 */
		int index2D;
		for (i = 1; i < r2; i++)
		{
			/* Process Row-i data 0 */
			index = k*r23 + i*r3;

			J = 0;
			if(i%blockEdgeSize==0)
				I++;
			realPrecision = pwrErrBound[K*R23+I*R3+J]; //J==0
			updateReqLength = 0;

			index2D = i*r3;
			pred2D = P0[index2D-r3] + P1[index2D] - P1[index2D-r3];
			diff = spaceFillingValue[index] - pred2D;

			itvNum = fabs(diff)/realPrecision + 1;

			if (itvNum < intvCapacity)
			{
				if (diff < 0) itvNum = -itvNum;
				type[index] = (int) (itvNum/2) + intvRadius;
				P0[index2D] = pred2D + 2 * (type[index] - intvRadius) * realPrecision;
			}
			else
			{
				if(updateReqLength==0)
				{
					computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
					reqBytesLength = reqLength/8;
					resiBitsLength = reqLength%8;
					updateReqLength = 1;
				}

				type[index] = 0;

				addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
				compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
				updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
				memcpy(preDataBytes,vce->curBytes,8);
				addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
				P0[index2D] = vce->data;
			}

			/* Process Row-i data 1 --> data r3-1 */
			for (j = 1; j < r3; j++)
			{
				index = k*r23 + i*r3 + j;

				if(j%blockEdgeSize==0)
				{
					J++;
					realPrecision = pwrErrBound[K*R23+I*R3+J];
					updateReqLength = 0;
				}
				index2D = i*r3 + j;
				pred3D = P0[index2D-1] + P0[index2D-r3]+ P1[index2D] - P0[index2D-r3-1] - P1[index2D-r3] - P1[index2D-1] + P1[index2D-r3-1];
				diff = spaceFillingValue[index] - pred3D;

				itvNum = fabs(diff)/realPrecision + 1;

				if (itvNum < intvCapacity)
				{
					if (diff < 0) itvNum = -itvNum;
					type[index] = (int) (itvNum/2) + intvRadius;
					P0[index2D] = pred3D + 2 * (type[index] - intvRadius) * realPrecision;
				}
				else
				{
					if(updateReqLength==0)
					{
						computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);
						reqBytesLength = reqLength/8;
						resiBitsLength = reqLength%8;
						updateReqLength = 1;
					}

					type[index] = 0;

					addDBA_Data(resiBitLengthArray, (unsigned char)resiBitsLength);
					compressSingleDoubleValue(vce, spaceFillingValue[index], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
					updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
					memcpy(preDataBytes,vce->curBytes,8);
					addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
					P0[index2D] = vce->data;
				}
			}
		}

		double *Pt;
		Pt = P1;
		P1 = P0;
		P0 = Pt;
	}
	if(r23!=1)
		free(P0);
	free(P1);
	int exactDataNum = exactLeadNumArray->size;

	TightDataPointStorageD* tdps;

	new_TightDataPointStorageD(&tdps, dataLength, exactDataNum,
			type, exactMidByteArray->array, exactMidByteArray->size,
			exactLeadNumArray->array,
			resiBitArray->array, resiBitArray->size,
			resiBitLengthArray->array, resiBitLengthArray->size,
			realPrecision, medianValue, (char)reqLength, quantization_intervals, pwrErrBoundBytes, pwrErrBoundBytes_size, radExpo);

	//free memory
	free_DBA(resiBitLengthArray);
	free_DIA(exactLeadNumArray);
	free_DIA(resiBitArray);
	free(type);

	convertTDPStoFlatBytes_double(tdps, newByteData, outSize);

	free(pwrErrBound);

	free(vce);
	free(lce);
	free_TightDataPointStorageD(tdps);
	free(exactMidByteArray);
}
