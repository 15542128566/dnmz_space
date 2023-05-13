#include<stdio.h>
#include<stdlib.h>

#define ABC(x) #x
#define DAY(x) myday##x

int main(){
    printf(ABC(Hello World!\n));
    int myday1 = 10;
    int myday2 = 20;
    // for (int i=1;i < 3; i++){
    //     printf("the day is %d\n", DAY(i));
    // }
    printf("the day is %d\n", DAY(1));
    printf("the day is %d\n", DAY(2));
    system("pause");

    return 0;
}