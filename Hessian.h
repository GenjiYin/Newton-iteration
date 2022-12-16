#include <iostream>
#include <windows.h>
#include <math.h>
using namespace std;

// 初始化函数, 传入向量的函数
// double func(double *arr, int len){
//     double out = 0;
//     for(int i=0; i<len; i++){
//         out += arr[i]*arr[i]*arr[i];
//     }
//     return out;
// }

// double *func_group(double *arr, int len){
//     // 方程传入数组
//     double x=arr[0];
//     double y=arr[1];
//     double z=arr[2];
//     double *p = new double[3];
//     p[0]=2*y*z+z*x-5*x*y-2;
//     p[1]=y*z-z*x+2*x*y-1;
//     p[2]=y*z-2*z*x+6*x*y-3;
//     return p;
// }

// 标量对向量求导
double *d(double f(double*, int), double *arr, int len){
    // arr为坐标
    double *p = new double[len];  // 用来放一阶导数的
    for(int i=0; i<len; i++){
        // 给当前分量加一个很小的值
        double a[len];
        for(int j=0; j<len; j++){
            if(j==i){
                a[j]=arr[j]+1e-6;
            }
            else{
                a[j]=arr[j];
            }
        }
        double f_1=f(a, len);
        double f_0=f(arr, len);
        p[i]=(f_1-f_0)/1e-6;
    }
    return p;
}

// 海塞矩阵
// 二阶导
double partial_sec(double f(double*, int), double *arr, int len, int index1, int index2, double epsilon_){
    // 对某两个维度求二阶偏导
    if(index1!=index2){
        double arr1[len];
        double arr2[len];
        double arr3[len];
        for(int i=0; i<len; i++){
            if(i==index1){
                arr1[i]=arr[i]+epsilon_;
                arr2[i]=arr[i];
                arr3[i]=arr[i]+epsilon_;
            }else if (i==index2)
            {
                arr1[i]=arr[i]+epsilon_;
                arr2[i]=arr[i]+epsilon_;
                arr3[i]=arr[i];
            }else{
                arr1[i]=arr[i];
                arr2[i]=arr[i];
                arr3[i]=arr[i];
            }
            
        }
        double out=f(arr1, len)-f(arr2, len)-f(arr3, len)+f(arr, len);
        return out/(epsilon_*epsilon_);
    }
    else{
        double arr1[len];
        double arr2[len];
        for(int i=0; i<len; i++){
            if(i==index1){
                arr1[i]=arr[i]+2*epsilon_;
                arr2[i]=arr[i]+epsilon_;
            }
            else{
                arr1[i]=arr[i];
                arr2[i]=arr[i];
            }
        }
        double out=f(arr1, len)-2*f(arr2, len)+f(arr, len);
        return out/(epsilon_*epsilon_);
    }
}

double **second(double f(double*, int), double *arr, int len){
    double **p=new double*[len];
    for(int i=0; i<len; i++){
        p[i]=new double[len];
    }
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            p[i][j]=partial_sec(f, arr, len, i, j, 1e-4);
        }
    }
    return p;
}

// 雅可比矩阵
// 向量函数对列向量变量求导
/*
df1/dx1 df2/dx1 ... dfn/dx1
...
df1/dxm df2/dfxm ...dfn/dfxm
*/
double **jacobi(double *func(double*, int), double *arr, int number, int len){
    // number: 方程组的数量 len: 坐标维度
    double **p = new double*[len];
    for(int i=0; i<len; i++){
        p[i]=new double[number];
    }
    for(int i=0; i<len; i++){
        for(int j=0; j<number; j++){
            double arr1[len];
            for(int o=0; o<len; o++){
                if(o==i){
                    arr1[o]=arr[o]+1e-6;
                }
                else{
                    arr1[o]=arr[o];
                }
            }
            double out1=func(arr, len)[j];
            double out2=func(arr1, len)[j];
            p[i][j]=(out2-out1)/1e-6;
        }
    }
    return p;
}
