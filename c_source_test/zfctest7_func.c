#include<stdio.h>
#include<stdlib.h>
#define SECONDS_OF_YEAR (356*24*3600UL)

int main(){
    printf("use printf func\n");

    int (*myshow)(const char *, ...);

    myshow = printf;
    myshow("use myshow as func printf\n");

    myshow("%d\n", SECONDS_OF_YEAR);

    system("pause");
    return 0;
}