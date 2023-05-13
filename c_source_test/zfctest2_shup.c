#include<stdio.h>
#include<stdlib.h>

#define ABC(x) #x

int main(){
    printf(ABC(Hello World!\n));
    system("pause");

    return 0;
}