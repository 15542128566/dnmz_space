#include<stdio.h>
#include<stdlib.h>

int main(){
    int a = 0x12345678;
    int b = 0x11223344;

    int *p1 = &a;
    int *p2 = &b;
    char *p3 = (char *)&b;

    printf("p1 is %x, p2 is %x, p2+1 is %x,%x\n", *p1, *p2,*(p2+1), p2[1]);
    printf("p3 is %x,p3+1 is %x\n", *p3, *(p3+1));

    system("pause");
}