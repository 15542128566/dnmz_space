#include<stdio.h>
#include<stdlib.h>

int main(){
    const int a = 0x12345678;
    int b = 0x11223344;

    int *p = &b;
    *(p+1) = 0x100;
    int *p2 =  &b;

    printf("a = %x,p2 =%x, p2+1=%x\n", a, *p2, *(p2+1));
    
    const int c = 0x12345678;
    int *p3 = &c;
    *p3 = 0x100;
    printf("c = %x\n", c);

    system("pause");

    return 0;
}