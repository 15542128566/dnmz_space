#include<stdio.h>
#include<stdlib.h>

// enum {MON=1,TUE,WED};
enum WEEK {MON=1,TUE,WED};

int main(){
    int a;
    printf("the size of a is %d\n",sizeof(a));
    printf("the size of a is %d\n",sizeof a);
    printf("the size of a is %lu\n",sizeof(a));
    printf("the size of a is %lu\n",sizeof a);

    printf("today is %d\n",MON);
    printf("today is %d\n",TUE);
    printf("today is %d\n",WED);

    enum WEEK week1;
    printf("week1 lenth is %lu\n",sizeof(week1));    

    int b = 10;
    b=~b;
    printf("b= %d\n",b);
    system("pause");

    return 0;
}