#include <math.h>
#include "fft.h"  

void c_abs(complex f[],double out[],int n)  
{  
  int i = 0;  
  double t;  
  for(i=0;i<n;i++)  
  {  
    t = f[i].real * f[i].real + f[i].imag * f[i].imag;  
    out[i] = sqrt(t);  
  }   
}    
  
void c_plus(complex a,complex b,complex *c)  
{  
  c->real = a.real + b.real;  
  c->imag = a.imag + b.imag;  
}  
  
void c_sub(complex a,complex b,complex *c)  
{  
  c->real = a.real - b.real;  
  c->imag = a.imag - b.imag;   
}  
  
void c_mul(complex a,complex b,complex *c)  
{  
  c->real = a.real * b.real - a.imag * b.imag;  
  c->imag = a.real * b.imag + a.imag * b.real;     
}  
  
void c_div(complex a,complex b,complex *c)  
{  
  c->real = (a.real * b.real + a.imag * b.imag)/(b.real * b.real +b.imag * b.imag);  
  c->imag = (a.imag * b.real - a.real * b.imag)/(b.real * b.real +b.imag * b.imag);  
} 


void cf_div(complex a,double b,complex *c)//增加复数与实数的除法
{
  c->real = a.real/b;
  c->imag = a.imag/b;
}

void cf_mul(complex a,double b,complex *c)//增加复数与实数的乘法
{
  c->real = a.real*b;
  c->imag = a.imag*b;
}

void Wn_i(int n,int i,complex *Wn,char flag)  
{  
  Wn->real = cos(2*PI*i/n);  
  if(flag == 1)  
  Wn->imag = -sin(2*PI*i/n);  
  else if(flag == 0)  
  Wn->imag = sin(2*PI*i/n);  
} 



//傅里叶变化  
void fft(int N,complex f[])  
{  
  complex t,wn={0,0};//中间变量  
  int i,j,k,m,n,l,r,M;  
  int la,lb,lc;  
  /*----计算分解的级数M=log2(N)----*/  
  for(i=N,M=1;(i=i/2)!=1;M++);   
  /*----按照倒位序重新排列原信号----*/  
  for(i=1,j=N/2;i<=N-2;i++)  
  {  
    if(i<j)  
    {  
      t=f[j];  
      f[j]=f[i];  
      f[i]=t;  
    }  
    k=N/2;  
    while(k<=j)  
    {  
      j=j-k;  
      k=k/2;  
    }  
    j=j+k;  
  }  
  
  /*----FFT算法----*/  
  for(m=1;m<=M;m++)  
  {  
    la=pow(2.0,m); //la=2^m代表第m级每个分组所含节点数       
    lb=la/2;    //lb代表第m级每个分组所含碟形单元数  
                 //同时它也表示每个碟形单元上下节点之间的距离  
    /*----碟形运算----*/  
    for(l=1;l<=lb;l++)  
    {  
      r=(l-1)*pow(2.0,M-m);     
      for(n=l-1;n<N-1;n=n+la) //遍历每个分组，分组总数为N/la  
      {  
        lc=n+lb;  //n,lc分别代表一个碟形单元的上、下节点编号       
        Wn_i(N,r,&wn,1);//wn=Wnr  
        c_mul(f[lc],wn,&t);//t = f[lc] * wn复数运算  
        c_sub(f[n],t,&(f[lc]));//f[lc] = f[n] - f[lc] * Wnr  
        c_plus(f[n],t,&(f[n]));//f[n] = f[n] + f[lc] * Wnr  
      }  
    }  
  }  
} 

void ifft(int N,complex f[])
{
	int i,j,k,m,n,l,r,M;
	int la,lb,lc;
	complex t={0,0};
	complex wn={0,0};

	/*----计算分解的级数M=log2(N)----*/
	for(i=N,M=1;(i=i/2)!=1;M++); 

	/*----将信号乘以1/N----*/
	for(i=0;i<N;i++) 
		cf_div(f[i],N,&f[i]);//f[i]=f[i]/complex(N,0.0);

	/*----IFFT算法----*/
	for(m=1;m<=M;m++)
	{
		la=pow(2.0,M+1-m); //la=2^m代表第m级每个分组所含节点数	
		lb=la/2;         //lb代表第m级每个分组所含碟形单元数
		                 //同时它也表示每个碟形单元上下节点之间的距离
		/*----碟形运算----*/
		for(l=1;l<=lb;l++)
		{
			r=(l-1)*pow(2.0,m-1);
			for(n=l-1;n<N-1;n=n+la) //遍历每个分组，分组总数为N/la
			{
				lc=n+lb;  //n,lc分别代表一个碟形单元的上、下节点编号     
				c_plus(f[n],f[lc],&t);      //t=f[n]+f[lc];
				Wn_i(N,r,&wn,0);
				c_sub(f[n],f[lc],&f[lc]);
				c_mul(f[lc],wn,&f[lc]);  //f[lc]=(f[n]-f[lc])*complex(cos(2*pi*r/N),sin(2*pi*r/N));
				f[n]=t;
			}
		}
	}

	/*----按照倒位序重新排列变换后信号----*/
	for(i=1,j=N/2;i<=N-2;i++)
	{
		if(i<j)
		{
			t=f[j];
			f[j]=f[i];
			f[i]=t;
		}
		k=N/2;
		while(k<=j)
		{
			j=j-k;
			k=k/2;
		}
		j=j+k;
	}
}