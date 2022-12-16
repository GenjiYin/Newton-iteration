#include "matrix.h"
#include "Hessian.h"

double *func_group(double *arr, int len){
    // 方程传入数组
    double x=arr[0];
    double y=arr[1];
    double z=arr[2];
    double *p = new double[3];
    p[0]=2*y*z+z*x-5*x*y-2;
    p[1]=y*z-z*x+2*x*y-1;
    p[2]=y*z-2*z*x+6*x*y-3;
    return p;
}

void *newtown(double *vector_func(double *, int), double *arr, int number, int dim, double factor, double epsilon_){
    /*
    参数数量和方程数量相同
    vector_func: 向量函数
    arr: 初始坐标
    number: 方程组的数量
    dim: 坐标维度
    factor: 阻尼因子
    epsilon: 阈值
    */
   double e=1000;
//    int c=0;
   while (e>epsilon_){
        // c++;
        double *F=func_group(arr, dim);
        // 行向量转化为列向量, 转化为二重指针
        double F_[dim][1];
        for(int i=0; i<dim; i++){
            F_[i][0]=F[i];
        }
        double *p[dim];
        for(int i=0; i<dim; i++){
            p[i]=F_[i];
        } 
        double **v = p;

        // 雅可比矩阵求逆
        double **J=jacobi(func_group, arr, number, dim);
        double **inv_J=inverse_matrix(J, dim);

        // 求下一轮迭代值
        double **temp=matrix_prod(inv_J, v, dim, number, 1);
        double delta[3];
        for(int i=0; i<3; i++){
            delta[i]=temp[i][0];
        }

        // 释放内存
        release_matrix(J, dim);
        release_matrix(inv_J, dim);
        for(int i=0; i<dim; i++){
            delete []temp[i];
        }
        delete []temp;

        // 数组相减
        for(int i=0; i<dim; i++){
            arr[i]=arr[i]-delta[i]*factor;
        }

        // 误差计算
        double *F__ = func_group(arr, dim);
        e=0;
        for(int i=0; i<dim; i++){
            e+=F__[i]*F__[i];
        }
        cout << e << endl;

        delete []F__;
        delete []F;
   }
}

int main(){
    double arr[]={1, 2, 3}; // 初始化值
    newtown(func_group, arr, 3, 3, 0.001, 2);    
    for(int i=0; i<3; i++){
        cout << arr[i]<<endl;
    }
    system("pause");
    return 0;
}