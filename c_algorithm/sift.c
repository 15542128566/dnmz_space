#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "./common.c"
#include "./bmp_util.c"
#include "./matrix_util.c"

#define MIN(a, b) ((a) < (b)) ? (a) : (b)
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define ABS(a) ((a) > 0.0) ? (a) : (-a)
#define SIFT_SIGMA 1.6          // 求sigma0的参数1
#define SIFT_INIT_SIGMA 0.5     // 假设的摄像头的尺度 默认取0.5
#define SIFT_SCALE_N 2          // 可求导最多差分高斯图片数
#define CONTRAST_THRESHOLD 0.04 // 阈值化参数T
#define EDGE_THRESHOLD 10.0     // 去除边缘效应阈值 gamma
#define SIFT_FIXPT_SCALE 255    // 像素范围
#define SIFT_MAX_INTERP_STEPS 5 // 寻找精确极值点时迭代最大次数
#define SIFT_IMG_BORDER 5       // 寻找精确极值点时查找像素点边界
#define PI 3.1415926
#define BIN_NUM 36 // 方向角分区

static const int DOG_SCALE_BASE = 3;

static int width, height; //原始图片的宽高

typedef struct point
{
    double posX;      // x坐标
    double posY;      // y坐标
    int posOS;        // 特征点所在的组信号与层序号
    double imgScale;  // 特征点所在组的图像尺寸
    double direction; // 特征点方向
} Point;

typedef struct point_list
{
    int length;
    Point *points;
} PointList;

void add_point(PointList *list, const Point *point)
{
    if (point == NULL)
    {
        return;
    }
    int len = list->length;
    if ((len % 10) == 0)
    {
        len += 10;
        list->points = (Point *)realloc(list->points, len * sizeof(Point));
    }
    list->points[list->length] = *point;
    list->length++;
}

// 求高斯金字塔
DoubleMat **get_guas_pyramid(unsigned char *img, int scaleN, double sigmaZero)
{
    int scale = scaleN + DOG_SCALE_BASE;
    int octave = (int)(log2(MIN(width, height)) - DOG_SCALE_BASE);
    double scaleK = pow(2, 1.0 / scaleN);
    // 高斯核矩阵 & 高斯金字塔
    double sigma[octave][scale];
    Mat *samplePyramid = (Mat *)malloc(octave * sizeof(Mat)); // 降采样图像数组
    // save_img_txt(img, width * height, 0, 0); // 存储原图像数据
    for (int o = 0; o < octave; ++o)
    {
        // 求高斯核矩阵
        for (int s = 0; s < scale; ++s)
        {
            sigma[o][s] = pow(scaleK, s) * sigmaZero * (1 << o);
        }
        // 降采样
        int step = (2 << o);
        samplePyramid[o] = create_mat(height / step, width / step);
        downsample(img, width, height, step, &samplePyramid[o]);
        save_mat_txt(&samplePyramid[o], o, 0); // 存中间过程数据：原图片采样数据
    }
    // 求高斯金字塔
    DoubleMat **guassianPyramid = (DoubleMat **)malloc(octave * scale * sizeof(DoubleMat *)); // 高斯图像数组
    int dim;
    for (int o = 0; o < octave; ++o)
    {
        guassianPyramid[o] = (DoubleMat *)malloc(scale * sizeof(DoubleMat));
        // 求高斯核矩阵
        for (int s = 0; s < scale; ++s)
        {
            dim = 2 * 3 * sigma[o][s] + 1; // 高斯核半径
            if ((dim % 2) == 0)
            {
                dim += 1; // 半径调整为奇数
            }
            // 高斯核矩阵
            DoubleMat *kernel = create_double_mat_ptr(dim, dim);
            guassian_kernel(sigma[o][s], dim, kernel);
            save_double_mat_txt(kernel, o * scale + s, 1); // 存中间过程数据：卷积核
            // 求卷积
            guassianPyramid[o][s] = create_double_mat(samplePyramid[o].rows, samplePyramid[o].cols);
            convolve(&samplePyramid[o], kernel, &guassianPyramid[o][s], 0, 0);
            save_double_mat_txt(&guassianPyramid[o][s], o * scale + s, 0); // 存中间过程数据：高斯金字塔数据
            release_double_mat_ptr(kernel);
            Mat *gausBmp = convert_double_mat_to_mat(&guassianPyramid[o][s]);
            save_mat_txt(gausBmp, o * scale + s, 1); // 存中间过程数据：高斯图金字塔像
            save_bmp_gray(o * scale + s, gausBmp->data, gausBmp->cols, gausBmp->rows);
            release_mat(gausBmp);
        }
    }

    return guassianPyramid;
}

// 求高斯差分金字塔
DoubleMat **get_dog(DoubleMat **guassianPyramid, int scaleN)
{
    int scale = scaleN + DOG_SCALE_BASE;
    int octave = (int)(log2(MIN(width, height)) - DOG_SCALE_BASE);
    // 求高斯差分金字塔
    DoubleMat **diffOfGuasPyramid = (DoubleMat **)malloc(octave * (scale - 1) * sizeof(DoubleMat *));
    for (int o = 0; o < octave; ++o)
    {
        diffOfGuasPyramid[o] = (DoubleMat *)malloc(scale * sizeof(DoubleMat));
        for (int s = 1; s < scale; ++s)
        {
            diffOfGuasPyramid[o][s - 1] = create_double_mat(guassianPyramid[o][s].rows, guassianPyramid[o][s].cols);
            for (int i = 0; i < guassianPyramid[o][s].rows; i++)
            {
                for (int j = 0; j < guassianPyramid[o][s].cols; j++)
                {
                    diffOfGuasPyramid[o][s - 1].data[i][j] = guassianPyramid[o][s].data[i][j] - guassianPyramid[o][s - 1].data[i][j];
                }
            }
            save_double_mat_txt(&diffOfGuasPyramid[o][s - 1], o * scale + s - 1, 2); // 存中间过程数据：高斯差分金字塔数据
        }
    }

    return diffOfGuasPyramid;
}

// 调整极值点使位置更精准
Point *adjust_local_extrema(DoubleMat **diffOfGuasPyramid, int o, int s, int x, int y, double sigma, int scaleN,
                            int *resX, int *resY, int *resS)
{
    double imgScale = 1.0 / SIFT_FIXPT_SCALE;
    double derivScale = imgScale * 0.5;       // 导数
    double secondDerivScale = imgScale;       // 二阶导
    double crossDerivScale = imgScale * 0.25; // 交叉求偏导

    DoubleMat img = diffOfGuasPyramid[o][s];
    DoubleMat imgPrev, imgNext;
    int iter = 0;
    double deriv[3];   // 一阶导数
    double xX, xY, xS; // 导数为零的解
    double dxx, dyy, dss, dxy, dxs, dys;
    while (iter < SIFT_MAX_INTERP_STEPS)
    {
        if (s < 1 || s > scaleN || x < SIFT_IMG_BORDER || (x >= img.rows - SIFT_IMG_BORDER) || y < SIFT_IMG_BORDER || (y >= img.cols - SIFT_IMG_BORDER))
        {
            return NULL;
        }
        imgPrev = diffOfGuasPyramid[o][s - 1];
        imgNext = diffOfGuasPyramid[o][s + 1];
        // 求导数
        deriv[0] = (img.data[x][y + 1] - img.data[x][y - 1]) * derivScale;
        deriv[1] = (img.data[x + 1][y] - img.data[x - 1][y]) * derivScale;
        deriv[2] = (imgNext.data[x][y] - imgPrev.data[x][y]) * derivScale;
        // 求二阶导数，构造hessian 矩阵
        double v2 = img.data[x][y] * 2;
        dxx = (img.data[x][y + 1] + img.data[x][y - 1] - v2) * secondDerivScale; // 对x方向求二阶偏导
        dyy = (img.data[x + 1][y] + img.data[x - 1][y] - v2) * secondDerivScale; // 对y方向求二阶偏导
        dss = (imgNext.data[x][y] + imgPrev.data[x][y] - v2) * secondDerivScale; // 对sigma方向求二阶偏导
        dxy = (img.data[x + 1][y + 1] - img.data[x + 1][y - 1] - img.data[x - 1][y + 1] + img.data[x - 1][y - 1]) * crossDerivScale;
        dxs = (imgNext.data[x][y + 1] - imgNext.data[x][y - 1] - imgPrev.data[x][y + 1] + imgPrev.data[x][y - 1]) * crossDerivScale;
        dys = (imgNext.data[x + 1][y] - imgNext.data[x - 1][y] - imgPrev.data[x + 1][y] + imgPrev.data[x - 1][y]) * crossDerivScale;
        double hessian[3][3] = {{dxx, dxy, dxs}, {dxy, dyy, dys}, {dxs, dys, dss}};
        // 求hessian逆矩阵
        double **hessianArr = creat_metrix_n3(hessian);
        double **hessianInverse = matrix_inver(hessianArr, 3);
        // 求出导数为零的解
        xX = -(hessianInverse[0][0] * deriv[0] + hessianInverse[0][1] * deriv[1] + hessianInverse[0][2] * deriv[2]);
        xY = -(hessianInverse[1][0] * deriv[0] + hessianInverse[1][1] * deriv[1] + hessianInverse[1][2] * deriv[2]);
        xS = -(hessianInverse[2][0] * deriv[0] + hessianInverse[2][1] * deriv[1] + hessianInverse[2][2] * deriv[2]);
        if ((ABS(xX) < 0.5) && ((ABS(xY) < 0.5)) && ((ABS(xS) < 0.5)))
        {
            break;
        }
        x += (int)round(xX);
        y += (int)round(xY);
        s += (int)round(xS);

        iter++;
    }
    if (iter >= SIFT_MAX_INTERP_STEPS)
    {
        return NULL;
    }
    if (s < 1 || s > scaleN || x < SIFT_IMG_BORDER || (x >= img.rows - SIFT_IMG_BORDER) || y < SIFT_IMG_BORDER || (y >= img.cols - SIFT_IMG_BORDER))
    {
        return NULL;
    }
    // 舍去对比度低的点
    double deltaVal = deriv[0] * xX + deriv[1] * xY + deriv[2] * xS;
    double resVal = img.data[x][y] * imgScale + deltaVal * 0.5;
    if ((ABS(resVal)) * scaleN < CONTRAST_THRESHOLD)
    {
        return NULL;
    }
    // 利用Hessian矩阵的迹和行列式计算主曲率的比值
    double tr = dxx + dyy;              // 求迹
    double det = dxx * dyy - dxy * dxy; // 求行列式
    if (det <= 0 || tr * tr * EDGE_THRESHOLD >= (EDGE_THRESHOLD + 1) * (EDGE_THRESHOLD + 1) * det)
    {
        return NULL;
    }
    Point *point = (Point *)malloc(sizeof(Point));
    point->posX = (x + xX) * (1 << o);
    point->posY = (y + xY) * (1 << o);
    point->posOS = o + (s << 8) + (((int)lround(xS + 0.5) * SIFT_FIXPT_SCALE)) << 16;
    point->imgScale = sigma * pow(2.0, (s + xS) / scaleN) * (1 << o) * 2;
    *resX = x;
    *resY = y;
    *resS = s;

    return point;
}

// 对梯度直方图做高斯平滑
void guas_smooth(double *hist, int binNum)
{
    double tmpHist[binNum + 4];
    for (int i = 0; i < binNum; i++)
    {
        tmpHist[i + 2] = hist[i];
    }
    tmpHist[0] = hist[binNum - 2];
    tmpHist[1] = hist[binNum - 1];
    tmpHist[binNum + 2] = hist[0];
    tmpHist[binNum + 3] = hist[1];
    for (int i = 0; i < binNum; i++)
    {
        hist[i] = (tmpHist[i] + tmpHist[i + 4]) * (1.0 / 16.0) +
                  (tmpHist[i + 1] + tmpHist[i + 3]) * (4.0 / 16.0) +
                  tmpHist[i + 2] * (6.0 / 16.0);
    }
}

// 确定关键点主方向
double *find_main_direction(DoubleMat img, int r, int c, int radius, double sigma, int binNum)
{
    double expfScale = -1.0 / (2.0 * sigma * sigma);
    double *hist = (double *)malloc(binNum * sizeof(double));
    for (int i = 0; i < binNum; i++)
    {
        hist[i] = 0.0;
    }

    int x, y, bin;
    double dx, dy, w, mag, angle;
    for (int i = -radius; i < radius + 1; i++)
    {
        y = r + i;
        if (y <= 0 || y >= img.rows - 1)
        {
            continue;
        }
        for (int j = -radius; j < radius + 1; j++)
        {
            x = c + j;
            if (x <= 0 || x >= img.cols - 1)
            {
                continue;
            }
            dx = img.data[y][x + 1] - img.data[y][x - 1];
            dy = img.data[y - 1][x] - img.data[y + 1][x];
            w = exp((i * i + j * j) * expfScale);
            // 计算梯度幅值大小
            mag = sqrt(dx * dx + dy * dy) / 1141 * SIFT_FIXPT_SCALE; // 1141 约为根号20，防止数据溢出
            // 计算梯度方向角
            angle = atan2(dy, dx) * 180 / PI;
            bin = (int)round((binNum / 360.0) * angle);
            if (bin >= binNum)
            {
                bin -= binNum;
            }
            if (bin < 0)
            {
                bin += binNum;
            }
            hist[bin] += (w * mag); // w为权重，离中心点越远权重越小
        }
    }
    // print_arr(hist, 36);
    // 高斯平滑
    guas_smooth(hist, binNum);

    return hist;
}

// 确定关键点位置
PointList *locate_key_points(DoubleMat **guassianPyramid, DoubleMat **diffOfGuasPyramid, double sigma, int scaleN)
{
    int scale = scaleN + DOG_SCALE_BASE;
    int octave = (int)(log2(MIN(width, height)) - DOG_SCALE_BASE);
    double oriSigma = 1.52;
    double oriRadius = 3 * oriSigma; // 找主方向搜寻范围
    double oriPeakRatio = 0.8;
    PointList *pointList;
    pointList->length = 0;
    pointList->points = (Point *)malloc(10 * sizeof(Point));
    Point *point;

    // 设定阈值，防止噪声影响
    double threshold = 0.5 * CONTRAST_THRESHOLD / (scaleN * SIFT_FIXPT_SCALE);
    DoubleMat imgPrev, img, imgNext; // 前一层、当前层、下一层 的高斯差分图像
    double val;                      // 当前中心点的值
    double neighbors[27] = {0};      // 前一层+当前层+下一层的相邻值， 共9+9+9个
    int valIsValid = 0;              // 判断是否满足阈值条件
    for (int o = 0; o < octave; ++o)
    {
        for (int s = 1; s < scale - 2; ++s)
        {
            imgPrev = diffOfGuasPyramid[o][s - 1];
            img = diffOfGuasPyramid[o][s];
            imgNext = diffOfGuasPyramid[o][s + 1];
            for (int i = 0; i < img.rows; i++)
            {
                for (int j = 0; j < img.cols; j++)
                {
                    // 先找出中心点的相邻值
                    int idx = 0;
                    for (int x = MAX(0, i - 1); x < (MIN(i + 2, img.rows)); x++)
                    {
                        for (int y = MAX(0, j - 1); y < (MIN(j + 2, img.cols)); y++)
                        {
                            neighbors[idx] = imgPrev.data[x][y];
                            neighbors[idx + 1] = img.data[x][y];
                            neighbors[idx + 2] = imgNext.data[x][y];
                            idx += 3;
                        }
                    }
                    // 判断是否满足阈值，且为极值点
                    val = img.data[i][j];
                    valIsValid = 0;
                    if ((ABS(val)) > threshold)
                    {
                        int isMax = 1, isMin = 1;
                        for (int k = 0; k < 27; k++)
                        {
                            if (val > 0 && val < neighbors[k])
                            {
                                isMax = 0;
                            }
                            else if (val < 0 && val > neighbors[k])
                            {
                                isMin = 0;
                            }
                        }
                        if ((val > 0 && isMax) || (val < 0 && isMin))
                        {
                            valIsValid = 1;
                        }
                    }
                    // 若满足阈值且为极值点，则寻找更为精确的极值点
                    if (valIsValid)
                    {
                        // printf("find extremum, o=%d, s=%d, i=%d, j=%d, val=%f\n", o, s, i, j, val);
                        int resX = i;
                        int resY = j;
                        int resS = s;
                        point = adjust_local_extrema(diffOfGuasPyramid, o, s, i, j, sigma, scaleN, &resX, &resY, &resS);
                        if (point == NULL)
                        {
                            continue;
                        }
                        // printf("before_x=%d,after_x=%d,before_y=%d,after_y=%d,before_s=%d,after_s=%d\n", i, resX, j, resY, s, resS);
                        double sclOctiv = point->imgScale * 0.5 / (1 << o);
                        double *hist = find_main_direction(guassianPyramid[o][resS], resX, resY, (int)round(oriRadius * sclOctiv), oriSigma * sclOctiv, BIN_NUM);
                        // print_arr(hist, 36);
                        double magMax = get_double_max(hist, BIN_NUM);
                        double magThreshold = magMax * oriPeakRatio;
                        for (int i = 0; i < BIN_NUM; i++)
                        {
                            int l = (i > 0) ? (i - 1) : (BIN_NUM - 1);
                            int r = (i < BIN_NUM - 1) ? (i + 1) : 0;
                            if (hist[i] > hist[l] && hist[i] > hist[r] && hist[i] >= magThreshold)
                            {
                                double bin = i + 0.5 * (hist[l] - hist[r]) / (hist[l] - 2.0 * hist[i] + hist[r]);
                                if (bin < 0)
                                {
                                    bin += BIN_NUM;
                                }
                                else if (bin >= BIN_NUM)
                                {
                                    bin -= BIN_NUM;
                                }
                                point->direction = (360.0 / BIN_NUM) * bin; // 计算特征点的方向
                                add_point(pointList, point);
                                printf("point: x=%f,y=%f,os=%d,scale=%f,direction=%f\n", point->posX, point->posY, point->posOS, point->imgScale, point->direction);
                            }
                        }
                    }
                }
            }
        }
    }

    return pointList;
}

// sift_n 取3
// sift_sigma 默认取1.6
// sift_init_sigma 假设的摄像头的尺度 默认取0.5
int sift(int sift_n, float sift_sigma, float sift_init_sigma)
{

    // 读图像文件的rpg
    unsigned char *imgRGB = read_rgb_bmp(&width, &height);
    if (imgRGB == NULL)
    {
        return 0;
    }
    // 转换灰度像素
    unsigned char *imgGray = rgbToGray(imgRGB, width, height);
    // 求高斯金字塔
    double sigmaZero = sqrt(sift_sigma * sift_sigma - sift_init_sigma * sift_init_sigma);
    DoubleMat **guasPyramid = get_guas_pyramid(imgGray, SIFT_SCALE_N, sigmaZero);
    // 求高斯差分金字塔
    DoubleMat **diffOfGuasPyramid = get_dog(guasPyramid, SIFT_SCALE_N);
    // 找极值点
    PointList *pointList = locate_key_points(guasPyramid, diffOfGuasPyramid, SIFT_SIGMA, SIFT_SCALE_N);

    printf("point_list_length=%d\n", pointList->length);
    // Mat *tmp = (Mat *)guassianPyramid[0];
    // save_bmp_gray("res.bmp", tmp->data, guassianPyramid[0][0].cols, guassianPyramid[0][0].rows);

    return 1;
}

int main()
{
    int res = sift(3, SIFT_SIGMA, SIFT_INIT_SIGMA);
    printf("sift res = %d\n", res);

    system("pause");

    return 0;
}