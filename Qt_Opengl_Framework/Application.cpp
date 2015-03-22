#include "Application.h"
#include "qt_opengl_framework.h"
//#include <vector>
//#include <iostream>
//#include <cstdlib>
//#include <cmath>

Application::Application()
{

}
Application::~Application()
{

}
void sort(int d[],int b,int t,int data[]){
	if(b<t){
		int n=b,i,temp;
        for(i=b;i<t;i++){
			if(d[i]>d[t]){
				temp=d[i];
				d[i]=d[n];
				d[n]=temp;
				temp=data[i];
				data[i]=data[n];
				data[n]=temp;
				n++;
			}
		}
		temp=d[t];
		d[t]=d[n];
		d[n]=temp;
		temp=data[i];
		data[i]=data[n];
		data[n]=temp;
        sort(d,b,n-1,data);
        sort(d,n+1,t,data);
	}
}
//****************************************************************************
//
// * 初始畫面，並顯示Ntust.png圖檔
// 
//============================================================================
void Application::createScene( void )
{
	
    ui_instance = Qt_Opengl_Framework::getInstance();
}

//****************************************************************************
//
// * 打開指定圖檔
// 
//============================================================================
void Application::openImage( QString filePath )
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * 刷新畫面
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
    ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * 畫面初始化
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * 儲存圖檔
// 
//============================================================================
void Application::saveImage(QString filePath )
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * 將圖檔資料轉換為RGB色彩資料
// 
//============================================================================
unsigned char* Application::To_RGB( void )
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (! img_data )
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0 ; j < img_width ; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j*4), rgb + (out_offset + j*3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB( unsigned char *rgba, unsigned char *rgb )
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0 ; i < 3 ; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i=0; i<img_height; i++)
	{
		for (int j=0; j<img_width; j++)
		{
			int offset_rgb = i*img_width*3+j*3;
			int offset_rgba = i*img_width*4+j*4;
			unsigned char gray = 0.3 * rgb[offset_rgb + rr] + 0.59 * rgb[offset_rgb + gg] + 0.11 * rgb[offset_rgb + bb];

			for (int k=0; k<3; k++)
				img_data[offset_rgba+k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();
	for (int i=0; i<img_height; i++)
	{
		for (int j=0; j<img_width; j++)
		{
			int offset_rgb = i*img_width*3+j*3;
			int offset_rgba = i*img_width*4+j*4;

            img_data[offset_rgba + rr] = rgb[offset_rgb + rr]&~31;//delete 5bit
            img_data[offset_rgba + gg] = rgb[offset_rgb + gg]&~31;//delete 5bit
            img_data[offset_rgba + bb] = rgb[offset_rgb + bb]&~64;//delete 6bit
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();
    //初始化統計用資料
    int histogram[32768]={};//記錄次數
    int histogramData[32768];//記錄顏色
	for(int i=0;i<32768;++i)
		histogramData[i]=i;
    //開始統計
	for (int i=0; i<img_height; i++){
		for (int j=0; j<img_width; j++){
			++histogram[
				((rgb[i*img_width*3+j*3 + rr]&~7)<<7) +
				((rgb[i*img_width*3+j*3 + gg]&~7)<<2) +
				((rgb[i*img_width*3+j*3 + bb]&~7)>>3)
            ];//去掉最後三個位元節省空間同時把三個顏色合成一個值
        }
    }
    //依數量多寡由大到小排
    sort(histogram,0,32767,histogramData);
    //決定每個pixel的顏色
	for (int i=0; i<img_height; i++){
		for (int j=0; j<img_width; j++){
			int offset_rgb = i*img_width*3+j*3;
            int offset_rgba = i*img_width*4+j*4;
            int nearColor=0;//最接近的顏色
            double len=262144;//跟最接近的顏色的距離
            for (int k=0;k<256;k++)	{//只與最高票的256色比
				int rl=((histogramData[k]>>7)&248)-rgb[offset_rgb + rr];
				int gl=((histogramData[k]>>2)&248)-rgb[offset_rgb + gg];
				int bl=((histogramData[k]<<3)&248)-rgb[offset_rgb + bb];
				double l=rl*rl+gl*gl+bl*bl;
                if(l<len){
					len=l;
					nearColor=histogramData[k];
				}
			}
            //把接成一個值的資料還原成三個值
			img_data[offset_rgba + rr] = ((nearColor>>7)&248);
			img_data[offset_rgba + gg] = ((nearColor>>2)&248);
            img_data[offset_rgba + bb] = ((nearColor<<3)&248);
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();
    for (int i=0; i<img_height*img_width; i++){
        //變灰
        unsigned char gray =
                0.30 * rgb[i*3 + rr] +
                0.59 * rgb[i*3 + gg] +
                0.11 * rgb[i*3 + bb];
        //變黑白
        img_data[i*4 + rr] = img_data[i*4 + gg] = img_data[i*4 + bb] = gray>127?255:0;
        img_data[i*4 + aa] = WHITE;
    }

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	unsigned char *rgb = this->To_RGB();
    for (int i=0; i<img_height*img_width; i++){
        //變灰
        unsigned char gray =
                0.30 * rgb[i*3 + rr] +
                0.59 * rgb[i*3 + gg] +
                0.11 * rgb[i*3 + bb];
        //取亂數(-50~50)
        int randNum = ((float)rand()/RAND_MAX-0.5)*100;
        //處理超界的值
        if(gray + randNum < 0)
            gray = 0;
        else if(gray + randNum > 255)
            gray = 255;
        else
            gray +=randNum;
        //變黑白
        img_data[i*4 + rr] = img_data[i*4 + gg] = img_data[i*4 + bb] = gray>127?255:0;
        img_data[i*4 + aa] = WHITE;
    }
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();
    double *gray = new double[img_height*img_width];//儲存灰階的值
    int direction = 1;//pixel前進的方向
    for (int i=0; i<img_height*img_width; i++)//先把灰階的值存進去
        gray[i] = 0.30 * rgb[i*3 + rr] + 0.59 * rgb[i*3 + gg] + 0.11 * rgb[i*3 + bb];
    for (int i=0; i<img_height; i++,direction*=-1){//到下一行就轉一次方向
        for (int j=(1-direction)/2*(img_width-1); j<img_width && j>=0; j+=direction){//依據方向算出起始點與前進
            int offset_gray = i*img_width+j;
            int offset_rgba = i*img_width*4+j*4;
            //算出顏色
            img_data[offset_rgba+rr] = gray[offset_gray]>127?255:0;
            //計算error
            double quantError = gray[offset_gray] - img_data[offset_rgba+rr];
            //把error分給其他pixel
            if(j+direction<img_width && j+direction>=0)//檢查邊界
                gray[offset_gray+direction] += quantError * 7 / 16 ;
            if(i!=img_height-1){
                gray[offset_gray+img_width] += quantError * 5 / 16 ;
                if(j+direction<img_width && j+direction>=0)
                    gray[offset_gray+direction+img_width] += quantError * 1 / 16 ;
                if(j-direction<img_width && j-direction>=0)
                    gray[offset_gray-direction+img_width] += quantError * 3 / 16 ;
            }
            //變黑白
            img_data[offset_rgba+gg] = img_data[offset_rgba+bb] = img_data[offset_rgba+rr];
            img_data[offset_rgba + aa] = WHITE;
        }
    }
    delete[] gray;
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();
    double whitePercent=0;//白色比率
    int histogram[256]={};//儲存每個顏色的數量
    //統計每個顏色的數量與白色的比例
	for (int i=0; i<img_height*img_width; i++){
            unsigned char gray =
                    0.30 * rgb[i*3 + rr] +
                    0.59 * rgb[i*3 + gg] +
                    0.11 * rgb[i*3 + bb];
			++histogram[gray];
            whitePercent+=gray;
    }
    whitePercent = whitePercent/(img_height*img_width)/255;
    int boundary = 255;//黑與白的邊界
    int pixelCount = 0;//數了幾個點
    //計算黑與白的邊界
    while(pixelCount<whitePercent*img_height*img_width)
        pixelCount+=histogram[boundary--];
    //變黑白
	for (int i=0; i<img_height*img_width; i++){
        unsigned char gray =
                0.30 * rgb[i*3 + rr] +
                0.59 * rgb[i*3 + gg] +
                0.11 * rgb[i*3 + bb];
        img_data[i*4 + rr] = img_data[i*4 + gg] = img_data[i*4 + bb] = gray>boundary?255:0;
        img_data[i*4 + aa] = WHITE;
    }

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
    unsigned char *rgb = this->To_RGB();
    int thresholdMatrix[16];//儲存Threshold Matrix值
    int randomMatrix[16];//決定儲存Threshold Matrix順序用
    //隨機取得Threshold Matrix順序
    for(int i=0;i<16;i++){//
        thresholdMatrix[i]=8+16*i;
        randomMatrix[i] = rand();
    }
    sort(randomMatrix, 0, 15, thresholdMatrix);
    //開始轉黑白
    for(int i=0;i<img_height;i+=4){
        for(int j=0;j<img_width;j+=4){
            for(int k=0;k<16;k++){
                //判斷是否有超界
                if(k%4+j<img_width&&k/4+i<img_height){
                    int offset = (i+k/4)*img_width+j+k%4;
                    img_data[offset*4+rr] =
                            img_data[offset*4+gg] =
                            img_data[offset*4+bb] =
                            0.30 * rgb[offset*3 + rr] +
                            0.59 * rgb[offset*3 + gg] +
                            0.11 * rgb[offset*3 + bb] >
                            thresholdMatrix[k]?255:0;
                    img_data[offset*4 + aa] = WHITE;
                }
            }
        }
    }

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();
    double *buffer = new double[img_height*img_width];//儲存顏色的值
    int direction = 1;//pixel前進的方向
    unsigned char mask[3]={192,224,224};//每個顏色要留下那些bit
    for(int colo=0;colo<3;colo++){//RGB分開做
        direction = 1;
        for (int i=0; i<img_height*img_width; i++)//先把顏色的值存進去
            buffer[i] = rgb[i*3 + colo];
        for (int i=0; i<img_height; i++,direction*=-1){//到下一行就轉一次方向
            for (int j=(1-direction)/2*(img_width-1); j<img_width && j>=0; j+=direction){//依據方向算出起始點與前進
                int offset = i*img_width+j;
                //算出顏色
                if(buffer[offset]>255||buffer[offset]<0)
                    img_data[offset*4+colo] = buffer[offset]>255?255:0;
                else
                    img_data[offset*4+colo] = ((unsigned char)buffer[offset])&mask[colo];
                //計算error
                double quantError = buffer[offset] - img_data[offset*4+colo];
                //把error分給其他pixel
                if(j+direction<img_width && j+direction>=0)//檢查邊界
                    buffer[offset+direction] += quantError * 7 / 16 ;
                if(i!=img_height-1){
                    buffer[offset+img_width] += quantError * 5 / 16 ;
                    if(j+direction<img_width && j+direction>=0)
                        buffer[offset+direction+img_width] += quantError * 1 / 16 ;
                    if(j-direction<img_width && j-direction>=0)
                        buffer[offset-direction+img_width] += quantError * 3 / 16 ;
                }
            }
        }
    }

    delete[] buffer;
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering( double filter[][5] )
{
	unsigned char *rgb = this->To_RGB();

	for(int x=0;x<img_height;++x){//pixel位置
		for(int y=0;y<img_width;++y){
			for(int color=0;color<3;++color){//顏色
				double sum = 0;
				for(int row=0;row<5;++row){//周圍位置
					for(int col=0;col<5;++col){
						int rowPosition = x-2+row;
						int colPosition = y-2+col;
						if(rowPosition<0)//檢查是否超界(使用反射)
							rowPosition = -rowPosition;
						if(rowPosition>=img_height)
							rowPosition = (img_height*2-2)-rowPosition;
						if(colPosition<0)
							colPosition = -colPosition;
						if(colPosition>=img_width)
							colPosition = (img_width*2-2)-colPosition;
						sum += rgb[(rowPosition*img_width+colPosition)*3+color]*filter[row][col];//把顏色乘以權重並加起來
					}
				}
				img_data[(x*img_width+y)*4+color] = sum;
			}
			img_data[(x*img_width+y)*4 + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::filtering( double **filter, int n )
{
	unsigned char *rgb = this->To_RGB();
//同上
	for(int x=0;x<img_height;++x){
		for(int y=0;y<img_width;++y){
			for(int color=0;color<3;++color){
				double sum = 0;
				for(int row=0;row<n;++row){
					for(int col=0;col<n;++col){
						int rowPosition = x-(n/2)+row;
						int colPosition = y-(n/2)+col;
						if(rowPosition<0)
							rowPosition = -rowPosition;
						if(rowPosition>=img_height)
							rowPosition = (img_height*2-(n/2))-rowPosition;
						if(colPosition<0)
							colPosition = -colPosition;
						if(colPosition>=img_width)
							colPosition = (img_width*2-(n/2))-colPosition;
						sum += rgb[(rowPosition*img_width+colPosition)*3+color]*filter[row][col];
					}
				}
				img_data[(x*img_width+y)*4+color] = sum;
				if(sum>255||sum<0)
					img_data[(x*img_width+y)*4+color]=sum>255?255:0;
			}
			img_data[(x*img_width+y)*4 + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	double filter[5][5] = {
		{1./25, 1./25, 1./25, 1./25, 1./25},
		{1./25, 1./25, 1./25, 1./25, 1./25},
		{1./25, 1./25, 1./25, 1./25, 1./25},
		{1./25, 1./25, 1./25, 1./25, 1./25},
		{1./25, 1./25, 1./25, 1./25, 1./25}
	};
    filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	double filter[5][5] = {
		{1/81., 2/81., 3/81., 2/81., 1/81.},
		{2/81., 4/81., 6/81., 4/81., 2/81.},
		{3/81., 6/81., 9/81., 6/81., 3/81.},
		{2/81., 4/81., 6/81., 4/81., 2/81.},
		{1/81., 2/81., 3/81., 2/81., 1/81.}
	};
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	double filter[5][5] = {
		{1/256., 4/256., 6/256., 4/256., 1/256.},
		{4/256., 16/256., 24/256., 16/256., 4/256.},
		{6/256., 24/256., 36/256., 24/256., 6/256.},
		{4/256., 16/256., 24/256., 16/256., 4/256.},
		{1/256., 4/256., 6/256., 4/256., 1/256.}
	};
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N( unsigned int n )
{
	if(n<1) return;//如果輸入小於1就不做
	//產生數列
	double *gaussianMask = new double[n];
	for(unsigned int i=0;i<n;++i){
		//算組合(Combination)
		gaussianMask[i]=1;
		for(unsigned int j=(n-1);j>i;--j)
			gaussianMask[i]*=j;
		for(unsigned int j=(n-1)-i;j>0;--j)
			gaussianMask[i]/=j;
	}
	//產生矩陣
	double** filter = new double*[n];
	for(unsigned int i=0;i<n;i++){
		filter[i] = new double[n];
		for(unsigned int j=0;j<n;j++)
			filter[i][j]=gaussianMask[i]*gaussianMask[j]/pow(pow((double)2,(int)(n-1)),2);
	}
	//使用濾鏡
	filtering(filter, n);
	//把矩陣與數列刪掉
	for(unsigned int i=0;i<n;i++)
		delete[] filter[i];
	delete[] filter;
	delete[] gaussianMask;
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{
	//轉灰階
	unsigned char *rgb = this->To_RGB();
	for (int i=0; i<img_height*img_width; i++)
		img_data[i*4 + rr] =
				img_data[i*4 + gg] =
				img_data[i*4 + bb] =
				0.30 * rgb[i*3 + rr] +
				0.59 * rgb[i*3 + gg] +
				0.11 * rgb[i*3 + bb];
	delete[] rgb;
	//用濾鏡
	double filter[5][5] = {
		{-1/256., -4/256., -6/256., -4/256., -1/256.},
		{-4/256., -16/256., -24/256., -16/256., -4/256.},
		{-6/256., -24/256., 220/256., -24/256., -6/256.},
		{-4/256., -16/256., -24/256., -16/256., -4/256.},
		{-1/256., -4/256., -6/256., -4/256., -1/256.}
	};
//	double filter[5][5] = {
//		{0, 0, 0, 0, 0},
//		{0, -1/16., -2/16., -1/16., 0},
//		{0, -2/16., 12/16., -2/16., 0},
//		{0, -1/16., -2/16., -1/16., 0},
//		{0, 0, 0, 0, 0},
//	};
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	double filter[5][5] = {
		{-1/256., -4/256., -6/256., -4/256., -1/256.},
		{-4/256., -16/256., -24/256., -16/256., -4/256.},
		{-6/256., -24/256., 476/256., -24/256., -6/256.},
		{-4/256., -16/256., -24/256., -16/256., -4/256.},
		{-1/256., -4/256., -6/256., -4/256., -1/256.}
	};
	filtering(filter);
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize( float scale )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate( float angleDegrees )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge( QString filePath )
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image( int tMethod )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::NPR_Paint_Layer( unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize )
{

}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke( const Stroke& s )
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) 
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) 
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height)) 
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared) 
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				} 
				else if (dist_squared == radius_squared + 1) 
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}





///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}
