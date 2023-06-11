#ifndef SIFT_H
#define SIFT_H
typedef struct mat
{
    int rows;
    int cols;
    unsigned char **data;
} Mat;

typedef struct double_mat
{
    int rows;
    int cols;
    double **data;
} DoubleMat;

typedef struct point
{
    double posX;      // x坐标
    double posY;      // y坐标
    int posOS;        // 特征点所在的组信号与层序号
    double imgScale;  // 特征点所在组的图像尺寸
    double direction; // 特征点方向
} Point;

typedef struct pointList
{
    int length;
    Point *point_list;
} PointList;

typedef struct descriptor
{
    int length;
    double *data;
} Descriptor;

typedef struct descriptorList
{
    int length;
    Descriptor *descriptor_list;
} DescriptorList;

#endif

void add_point(PointList *list, const Point *point);
void add_descriptor(DescriptorList *list, const Descriptor *descriptor);
Mat create_mat(int rows, int cols);
Mat *create_mat_ptr(int rows, int cols);
void release_mat(Mat *mat);
DoubleMat *create_double_mat_ptr(int rows, int cols);
DoubleMat create_double_mat(int rows, int cols);
void release_double_mat_ptr(DoubleMat *double_mat);
Mat **convert_double_mat_to_mat_ptr2(DoubleMat **src, int rows, int cols);
Mat *convert_double_mat_to_mat(const DoubleMat *src);
unsigned char *rgbToGray(unsigned char *src, int col, int row);
void downsample(unsigned char *img, int w, int h, int step, Mat *out);
void guassian_kernel(double sigma, int dim, DoubleMat *kernel);
void convolve(Mat *src, DoubleMat *kernel, DoubleMat *dest, int anchorx, int anchory);
unsigned char *read_rgb_bmp(int *width, int *height);
int custom_round(double number);
double get_double_max(const double *src, int length);

DoubleMat **get_guas_pyramid(unsigned char *img, int scaleN, double sigmaZero);
DoubleMat **get_dog(DoubleMat **guassianPyramid, int scaleN);
Point *adjust_local_extrema(DoubleMat **diffOfGuasPyramid, int o, int s, int x, int y, double sigma, int scaleN,
                            int *resX, int *resY, int *resS);
void guas_smooth(double *hist, int binNum);
double *find_main_direction(DoubleMat img, int r, int c, int radius, double sigma, int binNum);
PointList *locate_key_points(DoubleMat **guassianPyramid, DoubleMat **diffOfGuasPyramid, double sigma, int scaleN);
Descriptor *calcSIFTDescriptor(DoubleMat img, int x, int y, double angle, double scale);
DescriptorList *calcDescriptors(DoubleMat **guassianPyramid, PointList *pointList);
int sift(int sift_n, float sift_sigma, float sift_init_sigma);
