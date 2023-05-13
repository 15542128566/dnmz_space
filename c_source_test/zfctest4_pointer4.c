#include<stdio.h>
#include<stdlib.h>

int main (int argc,char **argv) {
    int i = 0;

    while (argv[i] != NULL) {
        printf("argv is %s\n", argv[i]);
        i++;
    };

    system("pause");

    return 0;
} 