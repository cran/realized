
#include "R.h"

#include "realized.h"




void justKernel(double *x, int *type, double *ans)
{
    ans[0] = K(x[0], type[0]); 

}
double K(double x, int type)
{
//    double thex = x[0];
    double thex = x;
	switch(type)
	{
      case 0:  // Rectangular
         return(1);
         break;
         
	  case 1: // Bartlett
	    return(1-thex);
	    break;
	   
	  case 2: //2nd Order 
	    return(1-2*thex+thex*thex);
	    break;
	    
	  case 3: //Epanechnikov
	    return(1-thex*thex);
	    break;

	  case 4: //Cubic
	    return(1-3*thex*thex+2*thex*thex*thex);
	    break;

	  case 5: //5th Order
	    return(1-10*thex*thex*thex + 15*thex*thex*thex*thex - 6*thex*thex*thex*thex*thex);
	    break;

	  case 6: //6th Order
	    return(1 - 15*thex*thex*thex*thex + 24*thex*thex*thex*thex*thex - 10*thex*thex*thex*thex*thex*thex);
	    break;

	  case 7: //7th Order
	    return(1 - 21*thex*thex*thex*thex*thex + 35*thex*thex*thex*thex*thex*thex - 15*thex*thex*thex*thex*thex*thex*thex);
	    break;

	  case 8: //8th Order
	    return(1 - 28*thex*thex*thex*thex*thex*thex + 48*thex*thex*thex*thex*thex*thex*thex - 21*thex*thex*thex*thex*thex*thex*thex*thex);
	    break;

	  case 9: //Parzen
	    if(thex > .5)
	    {
	        return( 2*(1-thex)*(1-thex)*(1-thex));
	    }
	    else
	    {
	    	return(1- 6*thex*thex +6*thex*thex*thex);
	    }
	    break;

	  case 10: //Tukey-Hanning
	    return((1 + sin(M_PI_2 - M_PI * thex))/2);
	    break;
	    
	  case 11: //Modified Tukey-Hanning
	    return((1 - sin(M_PI_2 - M_PI *(1 - thex) *(1 - thex)) )/2);
	    break;
	}
	return(-999);	
}


void kernel(double *a, int *na, int *q, int *adj, int *nw, double *ab, double *ans)
{
     	int i, j,lags = *q, nab = *na - 1;
        double theadj, nwadj;

        if(*nw == 1)
                lags = 2*lags;

        for(j = 0; j <= lags; j++)
		{        
	
        	for(i = 0; i <= nab-j; i++)
			     ab[j] += a[i] * a[i+j];
		}
        for(i = 0; i <= lags; i ++)
        {
              
              nwadj = 1;
              if(i > *q && *nw == 1)
              {
                  
                  nwadj = ((double)*q - ((double)i - (double)*q))/((double)*q); 
              }
               
              if(*adj == 0)
	     	 {
                  theadj = 1;
              }
              else
              {
                  theadj =  ((double)(nab+1)/((double)(nab+1)-(double)i));
	          }                                 
              if(i == 0)
	      	  {	
              		ans[0] += theadj * ab[i];        
              }
              else 
              {
	                ans[0] += 2 * nwadj *  theadj * ab[i];
              }
        }
}

void kernelEstimator(double *a, double *b, int *na, int *q, int *adj, int *type, double *ab,  double *ab2, double *ans)
{
     	int i, j,lags = *q, nab = *na - 1;
        double theadj, w;


        for(j = 0; j <= lags; j++)
		{        
	
        	for(i = 0; i <= nab-j; i++)
			     ab[j] += a[i] * b[i+j];
        	for(i = j; i <= nab; i++)
			     ab2[j] += a[i] * b[i-j];

		}
        for(i = 0; i <= lags; i ++)
        {
             if(i == 0)
             {
                 w = 1;
             }
             else
             { 
             	w = K( ((double)i-1.0)/(double)lags, type[0]);
             }	
              
              
             if(*adj == 0)
	     	 {
                  theadj = 1;
             }
             else
             {
                  theadj =  ((double)(nab+1)/((double)(nab+1)-(double)i));
	         }                                 
             
              if(i == 0)
	      	  {	
              		ans[0] += w * theadj * ab[i];        
              }
              else 
              {
	                ans[0] += w* ( theadj * ab[i] + theadj * ab2[i]);
              }
        }
}

void kernel2(double *a, double *b, int *na, int *q, int *adj, int *nw, double *ab,  double *ab2, double *ans)
{
     	int i, j,lags = *q, nab = *na - 1;
        double theadj, nwadj;

        if(*nw == 1)
                lags = 2*lags;

        for(j = 0; j <= lags; j++)
		{        
	
        	for(i = 0; i <= nab-j; i++)
			     ab[j] += a[i] * b[i+j];
        	for(i = j; i <= nab; i++)
			     ab2[j] += a[i] * b[i-j];

		}
        for(i = 0; i <= lags; i ++)
        {
              
              nwadj = 1;
              if(i > *q && *nw == 1)
              {
                  
                  nwadj = ((double)*q - ((double)i - (double)*q))/((double)*q); 
              }
               
              if(*adj == 0)
	     	 {
                  theadj = 1;
              }
              else
              {
                  theadj =  ((double)(nab+1)/((double)(nab+1)-(double)i));
	          }                                 
              if(i == 0)
	      	  {	
              		ans[0] += theadj * ab[i];        
              }
              else 
              {
	                ans[0] += nwadj *  theadj * ab[i];
	                ans[0] += nwadj *  theadj * ab2[i];
              }
        }
}



void subsample(double *a, double *b, int *na, int *m, int *period, double *tmpa, double *tmpb, int *tmpna, double *ans)
{
    int i,j,k;
    double starta, startb; 
    
    for(i = 0; i <= *period -1 ; i++)
    {
        starta = 0.0;
        startb = 0.0;
        for(j = 0; j < *tmpna; j++)
        {
             tmpa[j] = 0;
             tmpb[j] = 0;
        }
        for(j = 0; j <= *na-i-1; j++)
        {
        	tmpa[j/ *period] += a[j + i];
        	tmpb[j/ *period] += b[j + i];
        }
        for(j = 0; j < i; j++)
        {
            starta += a[j];
            startb += b[j];
        }
        for(k = 0; k < *tmpna; k++)
        {
            ans[i] += tmpa[k]*tmpb[k];
        }
        ans[i] += starta*startb;
	}
}

void rv(double *a, double *b, int *na, int *period, double *tmpa, double *tmpb, int *tmpna, double *ans)
{
    int j,k;
    
        for(j = 0; j <= *na-1; j++)
        {
        	tmpa[j/ *period] += a[j];
        	tmpb[j/ *period] += b[j];
        }
        for(k = 0; k < *tmpna; k++)
        {
            ans[0] += tmpa[k]*tmpb[k];
        }
}

void rvperiod(double *a, double *b, int *na, int *period, double *tmpa, double *tmpb, int *tmpna, double *ans)
{
    int j,k;
    
        for(j = 0; j <= *na-1; j++)
        {
        	tmpa[j/ *period] = a[j];
        	tmpb[j/ *period] = b[j];
        }
        for(k = 0; k < *tmpna; k++)
        {
            ans[0] += tmpa[k]*tmpb[k];
        }
}


void rfourth(double *a, double *b, int *na, int *period, double *tmpa, double *tmpb, int *tmpna, double *ans)
{
    int j,k;
    
        for(j = 0; j <= *na-1; j++)
        {
        	tmpa[j/ *period] += a[j];
        	tmpb[j/ *period] += b[j];
        }
        for(k = 0; k < *tmpna; k++)
        {
            ans[0] += tmpa[k]*tmpb[k]*tmpa[k]*tmpb[k];
        }
}

void rfourthlead(double *a, double *b, int *na, int *period, double *tmpa, double *tmpb, int *tmpna, double *ans)
{
    int j,k;
    
        for(j = 0; j <= *na-1; j++)
        {
        	tmpa[j/ *period] += a[j];
        	tmpb[j/ *period] += b[j];
        }
        for(k = 0; k < *tmpna-1; k++)
        {
            ans[0] += tmpa[k]*tmpb[k]*tmpa[k+1]*tmpb[k+1];
        }
}

void covkernel(double *a, double *b, int *na, int *q, int *adj, double *ans)
{
     	int i, j,lags = *q, nab = *na - 1;
        double theadj;

		for(i = lags; i <= nab-lags; i++)
        {
		    for(j =0; j<=lags; j++)
		       ans[0] += a[i] * b[i-j];
        }
        if(*adj == 0)
	    {
            theadj = 1;
        }
        else
        {
             theadj =  ((double)(nab+1)/((double)(nab+1)-(double)(2*lags)));
	    }                                 
        ans[0] = ans[0] * theadj;

}
	
void covcc(double *a, double *b, double *at, double *bt, int *na, int *nb, double *ans)
{
        int i, prevj = 0, j, nbi=*nb, nai = *na;
        double tmpRet=0;
        
        for(i = 0; i < nai; i++)
        {
            tmpRet=0;
            for(j=prevj; j < nbi; j++)
            {
                tmpRet += b[j];
                if(bt[j] > at[i])
                {
                    prevj = j;
                    j = nbi;
                }
                else if(bt[j] == at[i])
                {
                    prevj = j+1;
                    j = nbi;
                }
            }
            ans[i] = a[i]* tmpRet;
        } 
}
	
void pcovcc(double *a, double *ap, double *b, double *at,double *atp, double *bt, int *na, int *nap, int *nb, int *period, double *ans)
{
        int i, prevj = 0, j, nbi=*nb, nai = *na, napi = *nap;
        double tmpRet=0;
        
        for(i = 0; i < nai; i++)
        {
            ap[i/ *period] += a[i];
            atp[i/ *period] = at[i];
        }
        for(i = 0; i < napi; i++)
        {
            tmpRet=0;
            for(j=prevj; j < nbi; j++)
            {
                tmpRet += b[j];
                if(bt[j] > atp[i])
                {
                    prevj = j;
                    j = nbi;
                }
                else if(bt[j] == atp[i])
                {
                    prevj = j+1;
                    j = nbi;
                }
            }
            ans[i] = ap[i]* tmpRet;
        } 
}


void sametime(double *a, int *na, int *millis, double *tts)
{
    int prev = -1, i, j=0;
    for(i = 0; i < *na; i ++)
    {
    	if(millis[i] != prev)
    	{
    	    prev = millis[i];
    	    tts[j] = a[i];
    	    j = j +1;
    	}
    }
}	

void tocts(double *a, int *na, int *millis, int *millisstart, int *millisend, double *cts)
{
    int n = (*millisend - *millisstart)/1000 +1,i, j=0;
   for(i = 0; i < n; i ++)
    {
    	if((j < *na) && (millis[j] == ((i*1000) + *millisstart)))
    	{
    	    cts[i] = a[j];
    	    j = j+1;
    	}
    	else
    	{
    	   cts[i] = 0;
    	}
    }
 }
 
void portfolio(double *a, double *b, int *na, int *nb, int *millisa, int *millisb, int *millis,int *milliesn, double *port)
{
    int i, j=0, k=0;
    for(i = 0; i < *milliesn; i ++)
    {
    	if((j < *na) && (millisa[j] == millis[i]))
    	{
    	    port[i] = a[j];
    	    j = j+1;
    	}
    	else
    	{
    	   port[i] = 0;
    	}
        if((k < *nb) && (millisb[k] == millis[i]))
    	{
    	    port[i] += b[k];
    	    k = k+1;
    	}    	
    }
 }
 

