#include <stdio.h>
#include <stdlib.h>
// #define SWITCH

int main() {
    // 条件编译：如果定义了宏SWITCH则执行，否则不执行
#ifdef SWITCH
    printf("current msg: %s,%s,%d\n", __FUNCTION__, __FILE__, __LINE__);
#endif
    printf("Hello World!\n");
    system("pause");

    return 0;
}