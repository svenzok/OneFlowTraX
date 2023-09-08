/*
* find3x3peaks
*
* The calling syntax is:
*
*		peakList = find3x3peaks(imageStack,minIntensity,maxPeaks,ROIwidth)
*
* This is a MEX file for MATLAB.
*/
#include "mex.h"

/* declare main function routine */
void NMS(double *peakList, float *img,
         short nRows, short nCols, short nPages, 
         float minIntensity, unsigned long maxPeaks, short ROIwidth)
{
    short iCol, iRow, iPage, k, nWinner; //variable declaration
    unsigned long peakID = 0;
    int *canSkipNow;
    int *canSkipLater;
    int dCol[6] = {-1, 0, 1, -1,  0,  1};
    int dRow[6] = { 1, 1, 1, -1, -1, -1};

    canSkipNow = (int*)malloc(sizeof(int)*nCols);
    canSkipLater = (int*)malloc(sizeof(int)*nCols);
    for(iPage=0;iPage<nPages;++iPage)
    {
        for(iRow=(ROIwidth-1)/2; iRow<nRows-(ROIwidth-1)/2; ++iRow)
        {
            for(k=0; k<nCols-1; ++k)
            {
                /* make the skip info from the next line the current one
                and clear the next line skip info*/
                canSkipNow[k] = canSkipLater[k];
                canSkipLater[k] = 0;
            }

            for(iCol=(ROIwidth-1)/2; iCol<nCols-(ROIwidth-1)/2; ++iCol)
            {
                if ((canSkipNow[iCol])
                    || *(img + iRow*nCols + iCol + iPage*nCols*nRows) < minIntensity)
                {continue;} // skip pixel if masked
                if (  *(img + iRow*nCols + iCol + iPage*nCols*nRows)
                    > *(img + iRow*nCols + iCol + iPage*nCols*nRows + 1))
                {
                    canSkipNow[iCol + 1] = 1; // mask neighbor if current pixel is bigger
                }
                else {continue;}

                if (   *(img + iRow*nCols + iCol + iPage*nCols*nRows)
                    <= *(img + iRow*nCols + iCol + iPage*nCols*nRows - 1)) {continue;}

                /* if this line is reached, the current pixel is a line maximum*/
                nWinner = 0;
                for(k=0, nWinner=0; k<6; ++k, ++nWinner)
                {
                    /* visit neighborhood */
                    if (   *(img +  iRow           *nCols + iCol           + iPage*nCols*nRows)
                        <= *(img + (iRow + dRow[k])*nCols + iCol + dCol[k] + iPage*nCols*nRows))
                    {break;} // one of the neighborhood pixels is same or bigger

                }

                /* Mask future neighborhood pixels that were smaller than the
                current pixel. They are the first three in the
                neighborhood search, so the loop stops there. */
                for(k=0; k<3 && k<nWinner; ++k)
                {
                    canSkipLater[iCol + dCol[k]] = 1;
                }

                if (nWinner == 6)
                {
                    /* If true, then the pixel is a local intensity maximum
                    in a 3x3 neighborhood. Store the x and y pixel
                    coordinates indices pixel in MATLAB format. */
                    *(peakList+peakID) = iRow + 1; // x
                    *(peakList+peakID+maxPeaks) = iCol + 1; // y
                    *(peakList+peakID+2*maxPeaks) = iPage + 1; // frame
                    ++peakID; // next row in the output
                }
            }
        }
    }
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    float *img; //pointer to filtered image stack (single precision)
    short nRows,nCols,nPages; //size of image
    double *peakList; // pointer to the output array
    float minIntensity; // minimum pixel intensity to start processing
    unsigned long maxPeaks; // maximum number of expected peaks
    short ROIwidth; // width of ROI in pixels the peaks will be cut out later on
    const mwSize *dims; // to get the number of rows, columns, pages

    /*  check for proper number of arguments */
    if(nrhs!=4)
        mexErrMsgTxt("Four inputs required.");
    if(nlhs!=1)
        mexErrMsgTxt("One output required.");

    /*  create a pointer to the input matrix */
    img = mxGetData(prhs[0]);

    /*  get the dimensions of the matrix input */
    dims = mxGetDimensions(prhs[0]);
    nRows = dims[1]; // rows, y direction for an image
    nCols = dims[0]; // columns, x direction for an image
    nPages = dims[2]; // pages, frames for a image stack
    
    /* get the minimum intensity to process a pixel */
    minIntensity = mxGetScalar(prhs[1]);
        
    /* get the maximum number of peaks to process */
    maxPeaks = mxGetScalar(prhs[2]);

    /* get the ROI width */
    ROIwidth = mxGetScalar(prhs[3]);

    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateNumericMatrix(maxPeaks,3,mxDOUBLE_CLASS,mxREAL);

    /*  create a C pointer to a copy of the output matrix */
    peakList = mxGetDoubles(plhs[0]);

    /*  call the C subroutine */
    NMS(peakList,img,nRows,nCols,nPages,minIntensity,maxPeaks,ROIwidth);
}
