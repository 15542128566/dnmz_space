#include<stdio.h>
#include<stdlib.h>


struct abc {
    char a;
    short b;
    int c;
};

struct def {
    char d;
    int e;
    short f;
};

int main() {
    struct abc s1;
    struct def s2;

    printf("size of abc is %d\n", sizeof(s1));
    printf("size of def is %d\n", sizeof(s2));

    system("pause");

    return 0;
}