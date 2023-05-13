#include<stdio.h>
#include<stdlib.h>

int main() {
    char buf[12] = "abc";
    // strcpy(buf, "hello world");
    strncpy(buf, "hello world", sizeof(buf));

    printf("%s\n", buf);

    int a[10];
    int *p1 = a;
    printf("first address of a is %lu\n", p1);
    p1++;
    printf("first address of a is %lu\n", p1);

    int b[5][6];
    int (*p3)[6] = b;
    printf("first address of b is %lu\n", p3);
    printf("first address of b is %lu\n", &p3[0][1]);
    p3++;
    printf("first address of b is %lu\n", p3);

    system("pause");
    return 0;
}