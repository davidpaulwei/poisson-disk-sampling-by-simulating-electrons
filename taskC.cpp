#include "taskC.hpp"

//highlight image property
#define HIGHLIGHT_IMG_DEFAULT_COLOR Vec3b(0, 0, 0)
#define HIGHLIGHT_IMG_EDGE_COLOR Vec3b(255, 255, 255)
#define HIGHLIGHT_IMG_TAG_COLOR Scalar::all(0)

//tree image property
#define TREE_IMG_LENGTH 2048
#define TREE_IMG_HEIGHT 2048
#define TREE_IMG_BACKGROUND_COLOR Scalar::all(255)
#define TREE_IMG_NODE_DIA 64
#define TREE_IMG_NODE_COLOR Scalar::all(0)
#define TREE_IMG_LINE_COLOR Scalar::all(0)
#define TREE_IMG_LAYER_SPACING 32

//tree image layout functions
#define TREE_IMG_MAX_LAYER TREE_IMG_HEIGHT / (TREE_IMG_NODE_DIA + TREE_IMG_LAYER_SPACING)
#define TREE_IMG_LAYER_NODE_GAP(LAYER_NODE_COUNT) (float)(TREE_IMG_LENGTH - LAYER_NODE_COUNT * TREE_IMG_NODE_DIA) / LAYER_NODE_COUNT;
#define TREE_IMG_NODE_X(NODE_GAP, LAYER_NODE_INDEX) (NODE_GAP + TREE_IMG_NODE_DIA) * (LAYER_NODE_INDEX + 0.5);
#define TREE_IMG_NODE_Y(DEPTH) (TREE_IMG_NODE_DIA + TREE_IMG_LAYER_SPACING) * (DEPTH + 0.5);

extern int seedCount;

class GrfHeap {
public:
    GrfHeap(const Mat & mrkImage);
    ~GrfHeap();
    void sort();
    int findGrf(int minSize, int maxSize, int *& grfSizeList, int *& grfIdxList);
    int operator[](int i);
private:
    int binsearch(int grfSize) const;
    void rise(int ndIdx);
    void drop(int maxDropLen);
    inline void nodeSwap(int nodeIdxA, int nodeIdxB);
//data member:
    int * size;
    int * idx;
}; //Graph Heap Class

class BiTNode {
    friend BiTNode & huffmanCoding(int * const weightList, int len);
public:
    void showTree(Mat & treeImage);
    ~BiTNode();
private:
    BiTNode(int nodeWeight);
    BiTNode(BiTNode * const lchild, BiTNode * const rchild);
    int scanTree(int * nodeCount, int depth);
    void drawTree(Mat & treeImage, int * const nodeCount, int * nodeDrawn, int depth);
//data member:
    int weight;
    BiTNode * lchild;
    BiTNode * rchild;
};//Binary Tree Node Class

GrfHeap::GrfHeap(const Mat & mrkImage) : size(new int[seedCount]), idx(new int[seedCount]) {
    for(int i = 0; i < seedCount; size[i++] = 0);
    for(int i = 0; i < seedCount; idx[i++] = i);
    for(int i = 0; i < mrkImage.rows; i++) {
        for(int j = 0; j < mrkImage.cols; j++) {
            if(mrkImage.at<int>(i, j) > 0 && mrkImage.at<int>(i, j) <= seedCount)
                size[mrkImage.at<int>(i, j) - 1]++;
        }
    }
}

GrfHeap::~GrfHeap() {
    delete[] size;
    delete[] idx;
}

inline void GrfHeap::nodeSwap(int nodeIdxA, int nodeIdxB) {
    swap(size[nodeIdxA], size[nodeIdxB]);
    swap(idx[nodeIdxA], idx[nodeIdxB]);
}

void GrfHeap::rise(int ndIdx) {
    for(int i = ndIdx; i > 0 && size[i] > size[(i - 1) / 2]; nodeSwap(i, (i - 1) / 2), i = (i - 1) / 2);
}

void GrfHeap::drop(int maxDropLen) {
    int i = 0, maxChildIdx;
    while(true) {
        if(2 * i + 1 >= maxDropLen)
            break;
        else if(2 * i + 2 == maxDropLen)
            maxChildIdx = 2 * i + 1;
        else
            maxChildIdx = ((size[2 * i + 1] > size[2 * i + 2]) ? 2 * i + 1 : 2 * i + 2);
        if(size[i] >= size[maxChildIdx])
            break;
        nodeSwap(i, maxChildIdx);
        i = maxChildIdx;
    }
}

void GrfHeap::sort() {
    for(int i = 1; i < seedCount; i++)
        rise(i);
    for(int i = 1; i < seedCount; i++) {
        nodeSwap(0, seedCount - i);
        drop(seedCount - i);
    }
}

int GrfHeap::binsearch(int grfSize) const {
    int low = 0, high = seedCount;
    while(true) {
        if(high == low + 1)
            return low;
        if(size[(low + high) / 2] > grfSize)
            high = (low + high) / 2;
        else if(size[(low + high) / 2] == grfSize)
            return (low + high) / 2;
        else
            low = (low + high) / 2;
    }
}

int GrfHeap::findGrf(int minSize, int maxSize, int *& grfSizeList, int *& grfIdxList) {
    int minIdx = binsearch(minSize) + 1;
    int maxIdx = binsearch(maxSize);
    int grfCount = maxIdx - minIdx + 1;
    grfSizeList = new int[grfCount];
    grfIdxList = new int[grfCount];
    for(int i = minIdx; i <= maxIdx; i++) {
        grfSizeList[i - minIdx] = size[i];
        grfIdxList[i - minIdx] = idx[i];
    }
    return grfCount;
}

int GrfHeap::operator[](int i) {
    if(i >= 0 && i < seedCount)
        return size[i];
    else if(i == -1)
        return size[seedCount - 1];
    else
        return 0;
}

void highlight(Mat & hglImage, const Mat & marker, int * const grfIdxList, int * const grfTag, int grfCount) {
    Vec3b * colorTab = new Vec3b[seedCount];
    int * tagTab = new int[seedCount]();
    for(int i = 0; i < seedCount; i++)
        colorTab[i][0] = colorTab[i][1] = colorTab[i][2] = 127;
    for(int i = 0; i < grfCount; i++) {
        colorTab[grfIdxList[i]][0] = rand() % 255;
        colorTab[grfIdxList[i]][1] = rand() % 255;
        colorTab[grfIdxList[i]][2] = rand() % 255;
        tagTab[grfIdxList[i]] = grfTag[i];
    }
    for(int i = 0; i < marker.rows; i++) {
        for(int j = 0; j < marker.cols; j++) {
            int index = marker.at<int>(i, j);
            if(index == -1)
                hglImage.at<Vec3b>(i, j) = HIGHLIGHT_IMG_EDGE_COLOR;
            else if(index > 0 && index <= seedCount) {
                hglImage.at<Vec3b>(i, j) = colorTab[index - 1];
                if(tagTab[index - 1]) {
                    putText(hglImage, to_string(tagTab[index - 1]), Point(j, i), FONT_HERSHEY_PLAIN, 0.5, HIGHLIGHT_IMG_TAG_COLOR);
                    tagTab[index - 1] = 0;
                }
            } else
                hglImage.at<Vec3b>(i, j) = HIGHLIGHT_IMG_DEFAULT_COLOR;
        }
    }
    delete[] colorTab;
    delete[] tagTab;
}

int BiTNode::scanTree(int * nodeCount, int depth = 0) {
    nodeCount[depth]++;
    int lchildDepth, rchildDepth;
    if(lchild != NULL && depth < TREE_IMG_MAX_LAYER - 1)
        lchildDepth = lchild->scanTree(nodeCount, depth + 1);
    else
        lchildDepth = depth;
    if(rchild != NULL && depth < TREE_IMG_MAX_LAYER - 1)
        rchildDepth = rchild->scanTree(nodeCount, depth + 1);
    else
        rchildDepth = depth;
    if(lchildDepth > rchildDepth) {
        BiTNode * transChild = lchild;
        lchild = rchild;
        rchild = transChild;
    }
    return lchildDepth > rchildDepth ? lchildDepth : rchildDepth;
}

void BiTNode::drawTree(Mat & treeImage, int * const nodeCount, int * nodeDrawn, int depth = 0) {
    const float thisNodeGap = TREE_IMG_LAYER_NODE_GAP(nodeCount[depth]);
    const int thisNodeX = TREE_IMG_NODE_X(thisNodeGap, nodeDrawn[depth]);
    const int thisNodeY = TREE_IMG_NODE_Y(depth);
    nodeDrawn[depth]++;
    if(lchild != NULL && depth < TREE_IMG_MAX_LAYER - 1) {
        const float nextNodeGap = TREE_IMG_LAYER_NODE_GAP(nodeCount[depth + 1]);
        const int nextNodeX = TREE_IMG_NODE_X(nextNodeGap, nodeDrawn[depth + 1]);
        const int nextNodeY = TREE_IMG_NODE_Y(depth + 1);
        line(treeImage, Point(thisNodeX, thisNodeY), Point(nextNodeX, nextNodeY), TREE_IMG_LINE_COLOR);
        lchild->drawTree(treeImage, nodeCount, nodeDrawn, depth + 1);
    }
    if(rchild != NULL && depth < TREE_IMG_MAX_LAYER - 1) {
        const float nextNodeGap = TREE_IMG_LAYER_NODE_GAP(nodeCount[depth + 1]);
        const int nextNodeX = TREE_IMG_NODE_X(nextNodeGap, nodeDrawn[depth + 1]);
        const int nextNodeY = TREE_IMG_NODE_Y(depth + 1);
        line(treeImage, Point(thisNodeX, thisNodeY), Point(nextNodeX, nextNodeY), TREE_IMG_LINE_COLOR);
        rchild->drawTree(treeImage, nodeCount, nodeDrawn, depth + 1);
    }
    circle(treeImage, Point(thisNodeX, thisNodeY), TREE_IMG_NODE_DIA / 2, TREE_IMG_BACKGROUND_COLOR, -1);
    circle(treeImage, Point(thisNodeX, thisNodeY), TREE_IMG_NODE_DIA / 2, TREE_IMG_NODE_COLOR, 1);
    putText(treeImage, to_string(weight), Point(thisNodeX - 20, thisNodeY + 5), FONT_HERSHEY_PLAIN, 1, TREE_IMG_NODE_COLOR);
}

void BiTNode::showTree(Mat & treeImage) {
    int nodeCount[TREE_IMG_MAX_LAYER] = {0};
    int nodeDrawn[TREE_IMG_MAX_LAYER] = {0};
    scanTree(nodeCount);
    drawTree(treeImage, nodeCount, nodeDrawn);
}

BiTNode::BiTNode(int nodeWeight) : weight(nodeWeight), lchild(NULL), rchild(NULL) { }

BiTNode::BiTNode(BiTNode * lchildNode, BiTNode * rchildNode) : lchild(lchildNode), rchild(rchildNode), weight(lchildNode->weight + rchildNode->weight) { }

BiTNode::~BiTNode() {
    delete lchild;
    delete rchild;
    
}

BiTNode & huffmanCoding(int * weightList, int len) {
    BiTNode ** rootList = new BiTNode*[len];
    for(int i = 0; i < len; i++)
        rootList[i] = new BiTNode(weightList[i]);
    int minWeightRootIdxA, minWeightRootIdxB;
    int minWeightA, minWeightB;
    while(len > 1) {
        minWeightRootIdxA = 0;
        minWeightRootIdxB = 1;
        minWeightA = rootList[0]->weight;
        minWeightB = rootList[1]->weight;
        for(int i = 2; i < len; i++) {
            if(rootList[i]->weight > minWeightA && rootList[i]->weight > minWeightB)
                continue;
            if(minWeightA < minWeightB) {
                minWeightB = rootList[i]->weight;
                minWeightRootIdxB = i;
            } else {
                minWeightA = rootList[i]->weight;
                minWeightRootIdxA = i;
            }
        }
        rootList[minWeightRootIdxA] = new BiTNode(rootList[minWeightRootIdxA], rootList[minWeightRootIdxB]);
        rootList[minWeightRootIdxB] = rootList[len - 1];
        len--;
    }
    BiTNode & huffmanTreeRoot = **rootList;
    delete[] rootList;
    return huffmanTreeRoot;
}

int taskCMain(const Mat & marker) {
    //heap sort graphs by size
    GrfHeap grfSizeList = marker;
    double dTime = (double)getTickCount();
    grfSizeList.sort();
    dTime = (double)getTickCount() - dTime;
    cout << "[HEAP SORT] completed in " << dTime / getTickFrequency() << " seconds." << endl
         << "[HEAP SORT] " << seedCount << " graphs are now sorted by size,"
         << "with a minimum size of " << grfSizeList[0]
         << ", and a maximum size of " << grfSizeList[-1]
         << "." << endl;
    
    //highlight graphs whose size are within range
    int minSize, maxSize;
    cout << "[HIGHLIGHT GRAPH] Enter minimum and maximum size for highlighted graphs:" << endl << "- Minimum size: ";
    cin >> minSize;
    cout << "- Maximum size: ";
    cin >> maxSize;
    int * hglGrfSizeList = NULL, * hglGrfIdxList = NULL;
    int hglGrfCount = grfSizeList.findGrf(minSize, maxSize, hglGrfSizeList, hglGrfIdxList);
    if(hglGrfCount == 0) {
        cout << "[HIGHLIGHT GRAPH] No graphs found, end program." << endl;
        return 0;
    }
    Mat hglImage(marker.size(), CV_8UC3);
    highlight(hglImage, marker, hglGrfIdxList, hglGrfSizeList, hglGrfCount);
    cout << "[HIGHLIGHT GRAPH] " << hglGrfCount << " graphs found and highlighted." << endl;
    imshow("taskC", hglImage);
    waitKey(0);
    
    //huffman coding
    dTime = (double)getTickCount();
    BiTNode & huffmanTreeRoot = huffmanCoding(hglGrfSizeList, hglGrfCount);
    dTime = (double)getTickCount() - dTime;
    Mat huffmanTreeImage(TREE_IMG_HEIGHT, TREE_IMG_LENGTH, CV_8UC3, TREE_IMG_BACKGROUND_COLOR);
    huffmanTreeRoot.showTree(huffmanTreeImage);
    cout << "[HUFFMAN CODING] completed in " << dTime / getTickFrequency() << " seconds." << endl;
    imshow("taskC", huffmanTreeImage);
    waitKey(0);
    delete[] hglGrfSizeList;
    delete[] hglGrfIdxList;
    return 0;
}
