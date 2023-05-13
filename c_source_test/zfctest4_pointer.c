#include<stdio.h>
#include<stdlib.h>


int main(){
    int a = 0x12345678;
    char *p1 = &a;

    printf("value of pointer p is%x\n",*p1);

    // const
    char *p2 = "hello world!\n";
    printf("p2 first is %x\n", *p2);
    // *p2 = 'a';
    // printf("p2 first is %x\n", *p2);

    char buf[] = {"hello world\n"};
    char *p3 = buf;
    *p3 = 'a';
    *(p3 + 1)= 'b'; 
    printf("p3 first is %x\n", *p3);
    printf("p3 is %s\n", p3);
    // printf("p3 is %s\n", &p3);

    system("pause");
}