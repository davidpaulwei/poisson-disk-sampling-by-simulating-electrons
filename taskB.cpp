#include "taskB.hpp"

//hyperparams for 4-coloring algorithm
#define COLOR_COUNT 4
#define MAX_BACKTRACKING_LENGTH 5

//function params
#define SAME_COLOR_FOUND true
#define SAME_COLOR_NOT_FOUND false
#define ADJACENCY_ADDED true
#define ADJACENCY_EXISTED false
#define ADJACENCY_INVALID false
#define COLOR_UNSETTLED -1

//output image property
#define REGION_COLOR_1 Vec3b(255, 0, 0)
#define REGION_COLOR_2 Vec3b(255, 0, 255)
#define REGION_COLOR_3 Vec3b(0, 255, 0)
#define REGION_COLOR_4 Vec3b(0, 0, 255)
#define DEFAULT_COLOR Vec3b(0, 0, 0)
#define EDGE_COLOR Vec3b(255, 255, 255)
#define WATERSHED_IMG_TRANSPARENCY 0.5

extern int seedCount;

class GrfLNode {
    friend class GrfList;
public:
    GrfLNode();
    GrfLNode(int grfIndex, GrfLNode * const next);
    ~GrfLNode();
private:
    void add(int i);
//data member:
    int idx;
    GrfLNode * next;
}; //Adjacency Graph Link List Node Class


class GrfList {
public:
    GrfList(const Mat & marker);
    ~GrfList();
    int assignColor();
    void fillColor(const Mat & origImage, Mat & clredImage);
private:
    bool addAdjGrf(int i1, int i2);
    bool findSameColor(int index, int color);
//data member:
    GrfLNode * adjList; //Adjacency Graph Link List
    int * colorList; //Color List
}; //Graph List Class

GrfLNode::GrfLNode() : idx(-1), next(NULL) { }

GrfLNode::GrfLNode(int idx, GrfLNode * next) : idx(idx), next(next) { }

GrfLNode::~GrfLNode() {
    if(next != NULL)  delete next;
}

void GrfLNode::add(int i) {
    GrfLNode * newNode = new GrfLNode(i, next);
    next = newNode;
}

GrfList::GrfList(const Mat & marker) : adjList(new GrfLNode[seedCount]), colorList(new int[seedCount]) {
    for(int i = 0; i < seedCount; colorList[i++] = 0);
    for(int i = 1; i < marker.rows - 1; i++) {
        for(int j = 1; j < marker.cols - 1; j++) {
            if(marker.at<int>(i, j) != -1)
                continue;
            int firstGrfIdx = -1;
            if(marker.at<int>(i - 1, j) != -1)
                firstGrfIdx = marker.at<int>(i - 1, j);
            if(marker.at<int>(i + 1, j) != -1) {
                if(firstGrfIdx != -1 && marker.at<int>(i + 1, j) != firstGrfIdx) {
                    if(addAdjGrf(firstGrfIdx - 1, marker.at<int>(i + 1, j) - 1) == ADJACENCY_ADDED)
                        addAdjGrf(marker.at<int>(i + 1, j) - 1, firstGrfIdx - 1);
                    continue;
                }
                if(firstGrfIdx == -1)
                    firstGrfIdx = marker.at<int>(i + 1, j);
            }
            if(marker.at<int>(i, j - 1) != -1) {
                if(firstGrfIdx != -1 && marker.at<int>(i, j - 1) != firstGrfIdx) {
                    if(addAdjGrf(firstGrfIdx - 1, marker.at<int>(i, j - 1) - 1) == ADJACENCY_ADDED)
                        addAdjGrf(marker.at<int>(i, j - 1) - 1, firstGrfIdx - 1);
                    continue;
                }
                if(firstGrfIdx == -1)
                    firstGrfIdx = marker.at<int>(i, j - 1);
            }
            if(marker.at<int>(i, j + 1) != -1)
                if(firstGrfIdx != -1 && marker.at<int>(i, j + 1) != firstGrfIdx)
                    if(addAdjGrf(firstGrfIdx - 1, marker.at<int>(i, j + 1) - 1) == ADJACENCY_ADDED)
                        addAdjGrf(marker.at<int>(i, j + 1) - 1, firstGrfIdx - 1);
        }
    }
}

GrfList::~GrfList() {
    delete[] colorList;
    delete[] adjList;
}

bool GrfList::findSameColor(int grfIdx, int color) {
    GrfLNode * thisNode = adjList[grfIdx].next;
    while(true) {
        if(thisNode == NULL)
            return SAME_COLOR_NOT_FOUND;
        if(thisNode->idx > grfIdx)
            return SAME_COLOR_NOT_FOUND;
        if(colorList[thisNode->idx] == color)
            return SAME_COLOR_FOUND;
        thisNode = thisNode->next;
    }
}

bool GrfList::addAdjGrf(int thisGrfIdxIdx, int adjGrfIdx) {
    if(thisGrfIdxIdx < 0 || adjGrfIdx < 0)
        return ADJACENCY_INVALID;
    GrfLNode * thisNode = &(adjList[thisGrfIdxIdx]);
    while(true) {
        if(thisNode->next == NULL || (thisNode->next)->idx > adjGrfIdx) {
            thisNode->add(adjGrfIdx);
            return ADJACENCY_ADDED;
        }
        if((thisNode->next)->idx == adjGrfIdx)
            return ADJACENCY_EXISTED;
        thisNode = thisNode->next;
    }
}

int GrfList::assignColor() {
    int thisGrfIdx = 0, maxGrfIdx = 0, unstldGrfCount = 0;
    int colorHist[MAX_BACKTRACKING_LENGTH] = {0};
    while(thisGrfIdx < seedCount) {
        if(colorList[thisGrfIdx] == COLOR_COUNT) {
            colorList[thisGrfIdx] = 0;
            if(thisGrfIdx == maxGrfIdx - MAX_BACKTRACKING_LENGTH || colorList[thisGrfIdx - 1] == COLOR_UNSETTLED) {
                colorList[maxGrfIdx] = COLOR_UNSETTLED;
                for(int i = 0; i < MAX_BACKTRACKING_LENGTH; i++)
                    colorList[maxGrfIdx - MAX_BACKTRACKING_LENGTH + i] = colorHist[i];
                unstldGrfCount++;
                thisGrfIdx = ++maxGrfIdx;
            } else
                thisGrfIdx--;
            continue;
        }
        colorList[thisGrfIdx]++;
        if(findSameColor(thisGrfIdx, colorList[thisGrfIdx]) == SAME_COLOR_NOT_FOUND) {
            thisGrfIdx++;
            if(thisGrfIdx > maxGrfIdx) {
                maxGrfIdx = thisGrfIdx;
                for(int i = 0; i < MAX_BACKTRACKING_LENGTH; i++)
                    colorHist[i] = colorList[thisGrfIdx - MAX_BACKTRACKING_LENGTH + i];
            }
        }
    }
    for(thisGrfIdx = 0; thisGrfIdx < seedCount; thisGrfIdx++) {
        if(colorList[thisGrfIdx] == COLOR_UNSETTLED)
            colorList[thisGrfIdx] = rand() % COLOR_COUNT + 1;
    }
    return unstldGrfCount;
}

void GrfList::fillColor(const Mat & origImage, Mat & clredImage) {
    for(int i = 0; i < origImage.rows; i++) {
        for(int j = 0; j < origImage.cols; j++) {
            int index = origImage.at<int>(i, j);
            if(index == -1)
                clredImage.at<Vec3b>(i, j) = EDGE_COLOR;
            else if(index > 0 && index <= seedCount)
                switch(colorList[index - 1]) {
                    case 1: clredImage.at<Vec3b>(i, j) = REGION_COLOR_1; break;
                    case 2: clredImage.at<Vec3b>(i, j) = REGION_COLOR_2; break;
                    case 3: clredImage.at<Vec3b>(i, j) = REGION_COLOR_3; break;
                    case 4: clredImage.at<Vec3b>(i, j) = REGION_COLOR_4; break;
                    default: clredImage.at<Vec3b>(i, j) = DEFAULT_COLOR;
                }
            else
                clredImage.at<Vec3b>(i, j) = DEFAULT_COLOR;
        }
    }
}

void eraseRdnEdge(Mat & image) {
    for(int i = 1; i < image.rows - 1; i++) {
        for(int j = 1; j < image.cols - 1; j++) {
            if(image.at<Vec3b>(i, j) != EDGE_COLOR)
                continue;
            Vec3b firstGrfColor = EDGE_COLOR;
            if(image.at<Vec3b>(i - 1, j) != EDGE_COLOR)
                firstGrfColor = image.at<Vec3b>(i - 1, j);
            if(image.at<Vec3b>(i + 1, j) != EDGE_COLOR) {
                if(firstGrfColor == EDGE_COLOR)
                    firstGrfColor = image.at<Vec3b>(i + 1, j);
                else if(image.at<Vec3b>(i + 1, j) != firstGrfColor)
                    continue;
            }
            if(image.at<Vec3b>(i, j - 1) != EDGE_COLOR) {
                if(firstGrfColor == EDGE_COLOR)
                    firstGrfColor = image.at<Vec3b>(i, j - 1);
                else if(image.at<Vec3b>(i, j - 1) != firstGrfColor)
                    continue;
            }
            if(image.at<Vec3b>(i, j + 1) != EDGE_COLOR) {
                if(firstGrfColor == EDGE_COLOR)
                    firstGrfColor = image.at<Vec3b>(i, j + 1);
                else if(image.at<Vec3b>(i, j + 1) != firstGrfColor)
                    continue;
            }
            if(firstGrfColor != EDGE_COLOR)
                image.at<Vec3b>(i, j) = firstGrfColor;
        }
    }
}

int taskBMain(const Mat & marker) {
    //initialize graphs and graphlist
    GrfList grfGroup = marker;
    Mat srcImage = imread("Resources/lena.png", 1);
    Mat wsImage, grayImage;
    cvtColor(srcImage, wsImage, COLOR_BGR2GRAY);
    cvtColor(wsImage, grayImage, COLOR_GRAY2BGR);
    wsImage = Mat(srcImage.size(), CV_8UC3);
    
    //4-coloring
    double dTime = (double)getTickCount();
    grfGroup.assignColor();
    dTime = (double)getTickCount() - dTime;
    cout << "[4-COLORING] completed in " << dTime / getTickFrequency() << " seconds." << endl;
    
    //fill color & show image
    grfGroup.fillColor(marker, wsImage);
    eraseRdnEdge(wsImage);
    wsImage = grayImage * WATERSHED_IMG_TRANSPARENCY + wsImage * (1 - WATERSHED_IMG_TRANSPARENCY);
    imshow("taskB", wsImage);
    waitKey(0);
    return 0;
}
