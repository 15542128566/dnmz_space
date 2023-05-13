// #include <stdio.h>
#include "stdio.h"
#include<stdlib.h>

int isPrime(int p);

int main(){
    int a;
    int p;
    printf("请输入一个数字:\n");
    scanf("%d",&a);
    p = isPrime(a);
    printf("%d\n",p);
    system("pause");

    return 0;
}

int isPrime(int p){
    int flag;
    for (int i=2; i< p/2;i++){
        if (p%i == 0) {
            flag = 1;
            break;
        }
    }
    
    return flag;
}