#include <iostream>
#include <windows.h>
#include <math.h>
using namespace std;
#define epsilon 1e-15

/*
矩阵计算相关函数
inverse_matrix(arr, len): 矩阵求逆, 返回二重指针
输入二重指针与矩阵维度, 需释放内存

matrix_prod(arr1, arr2, row1, col1, col2): 矩阵相乘, 返回二重指针
输入两个二重指针 row1是第一个矩阵的行数, col1为第一个矩阵的列数也是第二个矩阵的行数, col2为第二个矩阵的列数 最后需释放内存

T(arr, row, col): 矩阵转置, 返回二重指针
输入二重指针 最后需释放内存

determinant(arr, len): 计算行列式, 返回浮点数
输入二重指针 不需要释放内存

release_matrix(arr, len): 释放二重指针内存
*/

// 数据结构
void release_matrix(double **arr, int len){
    // 释放堆区矩阵
    for(int i=0; i<len; i++){
        delete []arr[i];
    }
    delete []arr;
}

// 矩阵运算
// 行列式计算
double determinant(double **arr, int len){
    // 深拷贝
    double a[len][len]={0};
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            a[i][j]=arr[i][j];
        }
    }
    // 如果第一列全是0 返回零
    int row = 0;
    for(int i=0; i<len; i++){
        if(a[i][0]==0){
            row++;
        }
    }
    if(row==len-1){
        return 0;
    }
    
    // 化为上三角矩阵
    int number=0;
    for(int k=0; k<len-1; k++){
        int index = k;
        for(int i=k; i<len; i++){
            if(a[i][k]!=0){
                break;
            }
            else{
                index++;
            }
        }
        // if(index==len){abort();}
        if(index!=k){
            number++;
            for(int j=k; j<len; j++){
                double temp;
                temp=a[k][j];
                a[k][j]=a[index][j];
                a[index][j]=temp;
            }
        }
        for(int i=k+1; i<len; i++){
            double t=a[i][k];
            for(int j=k; j<len; j++){
                a[i][j]=a[i][j]-(t/a[k][k])*a[k][j]; // 系数要单列出来
            }
        }
    }
    double out =pow(-1, number);
    for(int i=0; i<len; i++){
        out *= a[i][i];
    }
    return out;
}

double **inverse_matrix(double **arr, int len){
    // 将矩阵作拷贝
    double **p=new double *[len];
    for(int i=0; i<len; i++){
        p[i]=new double[len];
    }
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            p[i][j]=arr[i][j];
        }
    }

    // 判断矩阵是否可逆
    double value=determinant(p, len);
    if(value<epsilon&&value>-epsilon){
        cout << "Uninversible" << endl;
        return 0;   //不可逆应该报错
    }
    release_matrix(p, len);
    
    // 初始化一个单位矩阵
    double **o=new double *[len];
    for(int i=0; i<len; i++){
        o[i]=new double[len];
    }
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            if(i==j){
                o[i][j]=1;
            }
            else{
                o[i][j]=0;
            }
        }
    }

    // 对原始矩阵进行深拷贝
    double a[len][len];
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            a[i][j]=arr[i][j];
        }
    }

    // 增广矩阵的方式对单位矩阵进行行变换
    // 上三角
    for(int k=0; k<len-1; k++){
        int index = k;
        for(int i=k; i<len; i++){
            if(a[i][k]!=0){
                break;
            }
            else{
                index++;
            }
        }
        // if(index==len){abort();}
        if(index!=k){
            for(int j=k; j<len; j++){
                // 原始矩阵交换位置
                double temp;
                temp=a[k][j];
                a[k][j]=a[index][j];
                a[index][j]=temp;

                // 单位矩阵交换位置
                double temp2;
                temp2=o[k][j];
                o[k][j]=o[index][j];
                o[index][j]=temp2;
            }
        }
        for(int i=k+1; i<len; i++){
            double t=a[i][k];
            for(int j=0; j<len; j++){
                a[i][j]=a[i][j]-(t/a[k][k])*a[k][j]; // 系数要单列出来
                o[i][j]=o[i][j]-(t/a[k][k])*o[k][j]; // 系数要单列出来
            }
        }
    }

    //下三角
    for(int i=len-2; i>=0; i--){
        double r[len-1-i];
        for(int j=i; j<len-1; j++){
            r[j-i]=(a[i][j+1]/a[j+1][j+1]);
        }
        for(int j=i+1; j<len; j++){ 
            for(int k=0; k<len; k++){
                o[i][k]=o[i][k]-r[j-i-1]*o[j][k];
            }
        }
    }

    // 每行归一
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            // a[i][j] = a[i][j]/a[i][i];
            o[i][j] = o[i][j]/a[i][i];
        }
    }
    return o;
}

// 矩阵乘法
double **matrix_prod(double **arr1, double **arr2, int row1, int col1, int col2){
    // 传入两个二维指针, row2==col1才能作矩阵相乘
    double **p=new double *[row1];
    for(int i=0; i<row1; i++){
        p[i]=new double[col2];
    }
    
    for(int i=0; i<row1; i++){
        for(int j=0; j<col2; j++){
            for(int k=0; k<col1; k++){
                p[i][j]+=arr1[i][k]*arr2[k][j];
            }
        }
    }
    return p;
}

void pr(double **arr, int len){
    for(int i=0; i<len; i++){
        for(int j=0; j<len; j++){
            cout << arr[i][j] <<" ";
        }
        cout << endl;
    }
}

// 矩阵转置
double **T(double **arr, int row, int col){
    double **p=new double*[col];
    for(int i=0; i<col; i++){
        p[i]=new double[row];
    }
    for(int i=0; i<col; i++){
        for(int j=0; j<row; j++){
            p[i][j]=arr[j][i];
        }
    }
    return p;
}