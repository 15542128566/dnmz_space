#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "./sift_util.c"
#include "./bmp_util.c"
#include "./matrix_util.c"
#include "sift.h"

#define MIN(a, b) ((a) < (b)) ? (a) : (b)
#define MAX(a, b) ((a) > (b)) ? (a) : (b)
#define ABS(a) ((a) > 0.0) ? (a) : (-a)
#define SIFT_SIGMA 1.6          // 求sigma0的参数1
#define SIFT_INIT_SIGMA 0.5     // 假设的摄像头的尺度 默认取0.5
#define SIFT_SCALE_N 3          // 可求导最多差分高斯图片数
#define CONTRAST_THRESHOLD 0.04 // 阈值化参数T
#define EDGE_THRESHOLD 10.0     // 去除边缘效应阈值 gamma
#define SIFT_FIXPT_SCALE 255    // 像素范围
#define SIFT_MAX_INTERP_STEPS 5 // 寻找精确极值点时迭代最大次数
#define SIFT_IMG_BORDER 5       // 寻找精确极值点时查找像素点边界
#define PI 3.1415926            // 圆周率
#define BIN_NUM 36              // 方向角分区
#define SIFT_DESCR_WIDTH 4      // 描述直方图的宽度
#define SIFT_DESCR_HIST_BINS 8  // 描述特征向量方向个数
#define SIFT_DESCR_SCL_FCTR 3.0 // 特征值范围系数
#define SIFT_INT_DESCR_FCTR 512.0

static const int DOG_SCALE_BASE = 3;
static int width, height; //原始图片的宽高
static int filter_cnt[7];

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
        int step = (1 << o);
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
    double deriv[3] = {0.0}; // 一阶导数
    double xX, xY, xS;       // 导数为零的解
    double dxx, dyy, dss, dxy, dxs, dys;
    while (iter < SIFT_MAX_INTERP_STEPS)
    {
        if (s < 1 || s > scaleN || x < SIFT_IMG_BORDER || (x >= img.rows - SIFT_IMG_BORDER) || y < SIFT_IMG_BORDER || (y >= img.cols - SIFT_IMG_BORDER))
        {
            filter_cnt[0]++;
            return NULL;
        }
        img = diffOfGuasPyramid[o][s];
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
        xY = -(hessianInverse[0][0] * deriv[0] + hessianInverse[0][1] * deriv[1] + hessianInverse[0][2] * deriv[2]);
        xX = -(hessianInverse[1][0] * deriv[0] + hessianInverse[1][1] * deriv[1] + hessianInverse[1][2] * deriv[2]);
        xS = -(hessianInverse[2][0] * deriv[0] + hessianInverse[2][1] * deriv[1] + hessianInverse[2][2] * deriv[2]);
        if ((ABS(xX) < 0.5) && ((ABS(xY) < 0.5)) && ((ABS(xS) < 0.5)))
        {
            filter_cnt[1]++;
            break;
        }
        x += (int)round(xX);
        y += (int)round(xY);
        s += (int)round(xS);

        iter++;
    }
    if (iter >= SIFT_MAX_INTERP_STEPS)
    {
        filter_cnt[2]++;
        *resX = x;
        *resY = y;
        *resS = s;
        return NULL;
    }
    if (s < 1 || s > scaleN || x < SIFT_IMG_BORDER || (x >= img.rows - SIFT_IMG_BORDER) || y < SIFT_IMG_BORDER || (y >= img.cols - SIFT_IMG_BORDER))
    {
        filter_cnt[3]++;
        return NULL;
    }
    // 舍去对比度低的点
    double deltaVal = deriv[0] * xY + deriv[1] * xX + deriv[2] * xS;
    double resVal = img.data[x][y] * imgScale + deltaVal * 0.5;
    if ((ABS(resVal)) * scaleN < CONTRAST_THRESHOLD)
    {
        filter_cnt[4]++;
        *resX = x;
        *resY = y;
        *resS = s;
        return NULL;
    }
    // 利用Hessian矩阵的迹和行列式计算主曲率的比值
    double tr = dxx + dyy;              // 求迹
    double det = dxx * dyy - dxy * dxy; // 求行列式
    if (det <= 0 || tr * tr * EDGE_THRESHOLD >= (EDGE_THRESHOLD + 1) * (EDGE_THRESHOLD + 1) * det)
    {
        filter_cnt[5]++;
        *resX = x;
        *resY = y;
        *resS = s;
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

    filter_cnt[6]++;
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
    PointList *pointList = (PointList *)malloc(sizeof(pointList));
    pointList->length = 0;
    pointList->point_list = (Point *)malloc(10 * sizeof(Point));
    Point *point;

    // 设定阈值，防止噪声影响
    double threshold = 0.5 * CONTRAST_THRESHOLD / (scaleN * SIFT_FIXPT_SCALE);
    DoubleMat imgPrev, img, imgNext; // 前一层、当前层、下一层 的高斯差分图像
    double val;                      // 当前中心点的值
    double neighbors[27] = {0};      // 前一层+当前层+下一层的相邻值， 共9+9+9个
    int valIsValid = 0;              // 判断是否满足阈值条件
    int cnt = 0;
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
                        cnt++;
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

    printf("cnt=%d\n", cnt);
    printf("filter_cnt_1 = %d\n", filter_cnt[0]);
    printf("filter_cnt_2 = %d\n", filter_cnt[1]);
    printf("filter_cnt_3 = %d\n", filter_cnt[2]);
    printf("filter_cnt_4 = %d\n", filter_cnt[3]);
    printf("filter_cnt_5 = %d\n", filter_cnt[4]);
    printf("filter_cnt_6 = %d\n", filter_cnt[5]);
    printf("filter_cnt_7 = %d\n", filter_cnt[6]);
    return pointList;
}

Descriptor *calcSIFTDescriptor(DoubleMat img, int x, int y, double angle, double scale)
{
    int d = SIFT_DESCR_WIDTH;
    int n = SIFT_DESCR_HIST_BINS;
    double cosT = cos(angle * (PI / 180));
    double sinT = sin(angle * (PI / 180));
    double bins_per_rad = n / 360.0;
    double exp_scale = -1.0 / (d * d * 0.5);
    double hist_width = SIFT_DESCR_SCL_FCTR * scale;
    double radius = (int)(round(hist_width * 1.414 * (d + 1) * 0.5)); // 求取值半径
    cosT /= hist_width;
    sinT /= hist_width;

    int max_size = (2 * radius + 1) * (2 * radius + 1);
    double *x_list = (double *)malloc(max_size * sizeof(double));
    double *y_list = (double *)malloc(max_size * sizeof(double));
    double *rbin_list = (double *)malloc(max_size * sizeof(double));
    double *cbin_list = (double *)malloc(max_size * sizeof(double));
    double *w_list = (double *)malloc(max_size * sizeof(double));
    double *mag_list = (double *)malloc(max_size * sizeof(double));
    double *angle_list = (double *)malloc(max_size * sizeof(double));
    int idx = 0;
    for (int i = -radius; i < radius + 1; i++)
    {
        for (int j = -radius; j < radius + 1; j++)
        {
            double c_rot = j * cosT - i * sinT;
            double r_rot = j * sinT - i * cosT;
            double rbin = r_rot + (int)(d / 2) - 0.5;
            double cbin = c_rot + (int)(d / 2) - 0.5;
            int r = y + i;
            int c = x + j;
            if (rbin > -1 && rbin < d && cbin > -1 && cbin < d && r > 0 && r < img.rows - 1 && c > 0 && c < img.cols - 1)
            {
                double dx = img.data[r][c + 1] - img.data[r][c - 1];
                double dy = img.data[r - 1][c] - img.data[r + 1][c];
                x_list[idx] = dx;
                y_list[idx] = dy;
                rbin_list[idx] = rbin;
                cbin_list[idx] = cbin;
                w_list[idx] = exp((c_rot * c_rot + r_rot * r_rot) * exp_scale);
                angle_list[idx] = atan2(dy, dx) * 180 / PI;
                mag_list[idx] = sqrt(dx * dx + dy * dy); // 1141 约为根号20，防止数据溢出;
                idx++;
            }
        }
    }
    save_descriptor_tmp_txt(x_list, idx, 0);
    save_descriptor_tmp_txt(y_list, idx, 1);
    save_descriptor_tmp_txt(rbin_list, idx, 2);
    save_descriptor_tmp_txt(cbin_list, idx, 3);
    save_descriptor_tmp_txt(w_list, idx, 4);
    save_descriptor_tmp_txt(mag_list, idx, 5);
    save_descriptor_tmp_txt(angle_list, idx, 6);
    if (idx == 0)
    {
        return NULL;
    }
    int len = idx;
    int hist_len = (d + 2) * (d + 2) * (n + 2);
    double *hist_list = (double *)malloc(hist_len * sizeof(double));
    for (int i = 0; i < hist_len; i++)
    {
        hist_list[i] = 0.0;
    }

    for (int k = 0; k < len; k++)
    {
        double rbin = rbin_list[k];
        double cbin = cbin_list[k];
        double obin = (angle_list[k] - angle) * bins_per_rad;
        double mag = mag_list[k] * w_list[k];
        int r0 = (int)rbin;
        int c0 = (int)cbin;
        int o0 = (int)obin;
        rbin -= r0;
        cbin -= c0;
        obin -= o0;
        if (o0 < 0)
        {
            o0 += n;
        }
        if (o0 >= n)
        {
            o0 -= n;
        }
        // 将直方图进行三线性插值更新
        double v_r1 = mag * rbin;
        double v_r0 = mag - v_r1;
        double v_rc11 = v_r1 * cbin;
        double v_rc10 = v_r1 - v_rc11;
        double v_rc01 = v_r0 * cbin;
        double v_rc00 = v_r0 - v_rc01;
        double v_rco111 = v_rc11 * obin;
        double v_rco110 = v_rc11 - v_rco111;
        double v_rco101 = v_rc10 * obin;
        double v_rco100 = v_rc10 - v_rco101;
        double v_rco011 = v_rc01 * obin;
        double v_rco010 = v_rc01 - v_rco011;
        double v_rco001 = v_rc00 * obin;
        double v_rco000 = v_rc00 - v_rco001;
        idx = ((r0 + 1) * (d + 2) + c0 + 1) * (n + 2) + o0;
        hist_list[idx] += v_rco000;
        hist_list[idx + 1] += v_rco001;
        hist_list[idx + (n + 2)] += v_rco010;
        hist_list[idx + (n + 3)] += v_rco011;
        hist_list[idx + (d + 2) * (n + 2)] += v_rco100;
        hist_list[idx + (d + 2) * (n + 2) + 1] += v_rco101;
        hist_list[idx + (d + 3) * (n + 2)] += v_rco110;
        hist_list[idx + (d + 3) * (n + 2) + 1] += v_rco111;
    }
    // 最终化直方图，方向直方图为环形
    Descriptor *dpr = (Descriptor *)malloc(sizeof(Descriptor));
    int dpr_len = d * d * n;
    dpr->length = dpr_len;
    dpr->data = (double *)malloc(dpr_len * sizeof(double));
    int iter = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            idx = ((i + 1) * (d + 2) + (j + 1)) * (n + 2);
            hist_list[idx] += hist_list[idx + n];
            hist_list[idx + 1] += hist_list[idx + n + 1];
            for (int k = 0; k < n; k++)
            {
                dpr->data[iter] = hist_list[idx + k];
                iter++;
            }
        }
    }
    save_descriptor_tmp_txt(hist_list, hist_len, 7);
    // for (int i = 0; i < hist_len; i++)
    // {
    //     printf("%f ", hist_list[i]);
    // }
    // printf("\n");

    // 将直方图复制到描述符中， 应用滞后阈值处理，并缩放结果，以便可以轻松地将其转换为字节数组
    double norm2 = 0.0;
    for (int k = 0; k < dpr_len; k++)
    {
        norm2 += dpr->data[k] * dpr->data[k];
    }
    double thr = sqrt(norm2) / 1141 * SIFT_FIXPT_SCALE;
    norm2 = 0.0;
    for (int k = 0; k < dpr_len; k++)
    {
        double val = MIN(dpr->data[k], thr);
        dpr->data[k] = val;
        norm2 += val * val;
    }
    norm2 = SIFT_INT_DESCR_FCTR / (MAX(sqrt(norm2), FLT_EPSILON));
    for (int k = 0; k < dpr_len; k++)
    {
        dpr->data[k] = MIN((MAX(dpr->data[k] * norm2, 0)), SIFT_FIXPT_SCALE);
    }

    return dpr;
}

DescriptorList *calcDescriptors(DoubleMat **guassianPyramid, PointList *pointList)
{
    int d = SIFT_DESCR_WIDTH;
    int n = SIFT_DESCR_HIST_BINS;
    int point_lens = pointList->length;
    PointList *new_point_list = (PointList *)malloc(sizeof(PointList));
    new_point_list->length = 0;
    new_point_list->point_list = (Point *)malloc(10 * sizeof(Point));
    DescriptorList *dpr_list = (DescriptorList *)malloc(sizeof(DescriptorList));
    dpr_list->length = 0;
    dpr_list->descriptor_list = (Descriptor *)malloc(10 * sizeof(Descriptor));
    Descriptor *dpr;
    for (int i = 0; i < point_lens; i++)
    {
        Point kpt = pointList->point_list[i];
        int o = kpt.posOS & SIFT_FIXPT_SCALE;        // 特征点所在组序号
        int s = (kpt.posOS >> 8) & SIFT_FIXPT_SCALE; // 特征点所在层序号
        double scale = 1.0 / (1 << o);               // 缩放倍数
        double size = kpt.imgScale * scale;          // 该特征点所在组的图像尺寸
        int x = (int)kpt.posX;
        int y = (int)kpt.posY;
        dpr = calcSIFTDescriptor(guassianPyramid[o][s], x, y, kpt.direction, scale);
        if (dpr != NULL)
        {
            add_descriptor(dpr_list, dpr);
            add_point(new_point_list, &kpt);
        }
    }
    *pointList = *new_point_list;

    return dpr_list;
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
    PointList *point_list = locate_key_points(guasPyramid, diffOfGuasPyramid, SIFT_SIGMA, SIFT_SCALE_N);
    printf("point_list_length=%d\n", point_list->length);
    DescriptorList *dpr_list = calcDescriptors(guasPyramid, point_list);
    printf("point_len=%d, descriptor_len=%d\n", point_list->length, dpr_list->length);
    save_key_points_txt(point_list, 0);
    save_descriptor_txt(dpr_list, 0);
    // 内存释放
    free(dpr_list);
    free(point_list);
    free(diffOfGuasPyramid);
    free(guasPyramid);
    free(imgGray);

    return 1;
}

int main()
{
    int res = sift(3, SIFT_SIGMA, SIFT_INIT_SIGMA);
    printf("sift res = %d\n", res);

    system("pause");

    return 0;
}