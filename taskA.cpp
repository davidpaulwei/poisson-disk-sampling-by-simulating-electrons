#include "taskA.hpp"

//physics engine hyperparams
#define IMG_LENGTH 512
#define IMG_WIDTH 512
#define CHARGE_UNIT_PER_ELEC_HIGH 0.05
#define CHARGE_UNIT_PER_ELEC_MED 0.03
#define CHARGE_UNIT_PER_ELEC_LOW 0.02
#define CHARGE_UNIT_PER_ELEC_XTRM_LOW 0.01
#define FRAME_WIDTH 1
#define MAX_REBOUND_DIST 2.
#define REBOUND_ACCURACY 10000.

//function params
#define MARK_AS_DOT 0
#define MARK_AS_CIRCLE 1

//watershed image property
#define DEFAULT_COLOR Vec3b(0, 0, 0)
#define EDGE_COLOR Vec3b(255, 255, 255)
#define DOT_COLOR Scalar(0, 255, 255)
#define WATERSHED_IMG_TRANSPARENCY 0.5

//basic property for watershed seed generation
int seedCount;
static int gridTtlCol;
static int gridTtlRow;
static float gridLength;
static float minDist;

class Seed {
    friend Seed * generateSeedGroup();
    friend void markSeedOnImage(Seed * seedGroup, Mat & marker, int flag);
public:
    Seed(Seed * fatherList);
    Seed();
    ~Seed();
private:
    void add(Seed * front);
    void remove();
    float distance(Point2f sd) const;
    void move(Seed * fatherList);
    void addForce(Seed & sd, float charge);
    void addForce(Point2f sd, float charge);
    int addTtlForce(Seed * fatherList, float charge);
//data member:
    Seed * prev;
    Point2f posVector;
    Point2f dispVector;
    Seed * next;
}; //Watershed Seed Double Link List Node Class

void fillRandomColor(const Mat & origImage, Mat & clredImage);

Seed::Seed(Seed * fatherList): posVector(rand() % IMG_LENGTH, rand() % IMG_WIDTH), prev(NULL), next(NULL) {
    add(fatherList + (int)((float)posVector.y / gridLength) * gridTtlCol + (int)((float)posVector.x / gridLength));
}

Seed::Seed(): prev(NULL), next(NULL) { }

Seed::~Seed() {
    if(next != NULL)  delete next;
}

void Seed::add(Seed * front) {
    next = front->next;
    if(next != NULL)
        next->prev = this;
    prev = front;
    front->next = this;
}

void Seed::remove() {
    prev->next = next;
    if(next != NULL)
        next->prev = prev;
}

float Seed::distance(Point2f sd) const {
    return sqrt(pow(posVector.x - sd.x, 2) + pow(posVector.y - sd.y, 2));
}

void Seed::move(Seed * fatherList) {
    if(dispVector == Point2f(0, 0))
        return;
    int prevGridIndex = (int)((float)posVector.y / gridLength) * gridTtlCol + (int)((float)posVector.x / gridLength);
    posVector += dispVector;
    dispVector = Point2f(0, 0);
    if(posVector.x > IMG_LENGTH - 1 - FRAME_WIDTH)
        posVector.x = IMG_LENGTH - 1 - FRAME_WIDTH - MAX_REBOUND_DIST * (rand() % (int)REBOUND_ACCURACY) / REBOUND_ACCURACY;
    else if (posVector.x < FRAME_WIDTH)
        posVector.x = FRAME_WIDTH + MAX_REBOUND_DIST * (rand() % (int)REBOUND_ACCURACY) / REBOUND_ACCURACY;
    if(posVector.y > IMG_WIDTH - 1 - FRAME_WIDTH)
        posVector.y = IMG_WIDTH - 1 - FRAME_WIDTH - MAX_REBOUND_DIST * (rand() % (int)REBOUND_ACCURACY) / REBOUND_ACCURACY;
    else if (posVector.y < FRAME_WIDTH)
        posVector.y = FRAME_WIDTH + MAX_REBOUND_DIST * (rand() % (int)REBOUND_ACCURACY) / REBOUND_ACCURACY;
    int gridIndex = (int)((float)posVector.y / gridLength) * gridTtlCol + (int)((float)posVector.x / gridLength);
    if (prevGridIndex == gridIndex)
        return;
    remove();
    add(fatherList + gridIndex);
}

void Seed::addForce(Seed & sd, float charge) {
    const float dist = distance(sd.posVector);
    if(dist == 0)
        return;
    if(dist >= minDist)
        return;
    const float dispXAxis = (posVector.x - sd.posVector.x) * charge / pow(dist, 3);
    const float dispYAxis = (posVector.y - sd.posVector.y) * charge / pow(dist, 3);
    dispVector.x += dispXAxis;
    dispVector.y += dispYAxis;
    sd.dispVector.x -= dispXAxis;
    sd.dispVector.y -= dispYAxis;
}

inline void Seed::addForce(Point2f sd, float charge) {
    const float dist = distance(sd);
    if(dist >= minDist)
        return;
    dispVector.x += (posVector.x - sd.x) * charge / pow(dist, 3);
    dispVector.y += (posVector.y - sd.y) * charge / pow(dist, 3);
}

int Seed::addTtlForce(Seed * fatherList, float charge) {
    const int col = (int)((float)posVector.x / gridLength), row = (int)((float)posVector.y / gridLength);
    const int gridIndex = row * gridTtlCol + col;
    Seed * cmpList = this;
    int unstldFlag = 0;
    while(cmpList->next != NULL) {
        cmpList = cmpList->next;
        addForce(*cmpList, charge);
    }
    //add real force (only from points behind)
    if(col != gridTtlCol - 1 && fatherList[gridIndex + 1].next != NULL) {
        cmpList = fatherList[gridIndex + 1].next;
        while (true) {
            addForce(*cmpList, charge);
            if(cmpList->next == NULL)
                break;
            cmpList = cmpList->next;
        }
    }
    if(row != gridTtlRow - 1) {
        if(col != 0 && fatherList[gridIndex + gridTtlCol - 1].next != NULL) {
            cmpList = fatherList[gridIndex + gridTtlCol - 1].next;
            while (true) {
                addForce(*cmpList, charge);
                if(cmpList->next == NULL)
                    break;
                cmpList = cmpList->next;
            }
        }
        if(fatherList[gridIndex + gridTtlCol].next != NULL) {
            cmpList = fatherList[gridIndex + gridTtlCol].next;
            while (true) {
                addForce(*cmpList, charge);
                if(cmpList->next == NULL)
                    break;
                cmpList = cmpList->next;
            }
        }
        if(col != gridTtlCol - 1 && fatherList[gridIndex + gridTtlCol + 1].next != NULL) {
            cmpList = fatherList[gridIndex + gridTtlCol + 1].next;
            while (true) {
                addForce(*cmpList, charge);
                if(cmpList->next == NULL)
                    break;
                cmpList = cmpList->next;
            }
        }
    }
    if(dispVector != Point2f(0, 0))
        unstldFlag = 1;
    //add shadow point force
    if(row == 0) {
        addForce(Point2f(col * gridLength, -gridLength / 2 * sqrt(3)), charge);
        addForce(Point2f((col + 1) * gridLength, -gridLength / 2 * sqrt(3)), charge);
    }
    if(row == gridTtlRow) {
        addForce(Point2f(col * gridLength, IMG_WIDTH + gridLength / 2 * sqrt(3)), charge);
        addForce(Point2f((col + 1) * gridLength, IMG_WIDTH + gridLength / 2 * sqrt(3)), charge);
    }
    if(col == 0) {
        addForce(Point2f(-gridLength / 2 * sqrt(3), row * gridLength), charge);
        addForce(Point2f(-gridLength / 2 * sqrt(3), (row + 1) * gridLength), charge);
    }
    if(col == gridTtlCol) {
        addForce(Point2f(IMG_LENGTH + gridLength / 2 * sqrt(3), row * gridLength), charge);
        addForce(Point2f(IMG_LENGTH + gridLength / 2 * sqrt(3), (row + 1) * gridLength), charge);
    }
    return unstldFlag;
}

Seed * generateSeedGroup() {
    Seed * fatherList = new Seed[gridTtlCol * gridTtlRow]();
    const float minDist = sqrt(IMG_WIDTH * IMG_LENGTH / (float)seedCount);
    const float unitCharge = pow(minDist, 3);
    float charge = unitCharge * CHARGE_UNIT_PER_ELEC_HIGH;
    //initialize random seeds.
    for(int i = 0; i < seedCount; i++)
        new Seed(fatherList);
    //move seeds until fit requirements.
    int unstldSeedCount;
    int iterCount = 0;
    Seed * thisList;
    do {
        unstldSeedCount = 0;
        //add force
        for(int i = 0; i < gridTtlCol * gridTtlRow; i++) {
            if (fatherList[i].next == NULL)
                continue;
            thisList = fatherList[i].next;
            while (true) {
                unstldSeedCount += thisList->addTtlForce(fatherList, charge);
                if(thisList->next == NULL)
                    break;
                thisList = thisList->next;
            }
        }
        //move
        for(int i = 0; i < gridTtlCol * gridTtlRow; i++) {
            if (fatherList[i].next == NULL)
                continue;
            thisList = fatherList[i].next;
            while (true) {
                    thisList->move(fatherList);
                if(thisList->next == NULL)
                    break;
                thisList = thisList->next;
            }
        }
        //change charge const
        if(unstldSeedCount < seedCount / 20)
            charge = unitCharge * CHARGE_UNIT_PER_ELEC_XTRM_LOW;
        else if (unstldSeedCount < seedCount / 3)
            charge = unitCharge * CHARGE_UNIT_PER_ELEC_LOW;
        else if (unstldSeedCount < seedCount * 2 / 3)
            charge = unitCharge * CHARGE_UNIT_PER_ELEC_MED;
        else
            charge = unitCharge * CHARGE_UNIT_PER_ELEC_HIGH;
        iterCount++;
    } while(unstldSeedCount);
    cout << "[PHYSICS ENGINE] Electrons stabilized at iteration count " << iterCount << "." << endl;
    return fatherList;
}

void markSeedOnImage(Seed * seedGroup, Mat & marker, int flag = MARK_AS_DOT){
    Seed * thisList;
    int pointIdx = 0;
    for(int i = 0; i < gridTtlCol * gridTtlRow; i++) {
        if (seedGroup[i].next == NULL)
            continue;
        thisList = seedGroup[i].next;
        while (true) {
            if(flag == MARK_AS_DOT)
                marker.at<int>(thisList ->posVector) = ++pointIdx;
            else if(flag == MARK_AS_CIRCLE)
                circle(marker, Point(thisList->posVector.x, thisList->posVector.y), 1, DOT_COLOR, -1);
            if(thisList->next == NULL)
                break;
            thisList = thisList->next;
        }
    }
}

void fillRandomColor(const Mat & origImage, Mat & clredImage) {
    Vec3b * colorTab = new Vec3b[seedCount];
    for(int i = 0; i < seedCount; i++) {
        colorTab[i][0] = rand() % 255;
        colorTab[i][1] = rand() % 255;
        colorTab[i][2] = rand() % 255;
    }
    for(int i = 0; i < origImage.rows; i++) {
        for(int j = 0; j < origImage.cols; j++) {
            int index = origImage.at<int>(i, j);
            if(index == -1)
                clredImage.at<Vec3b>(i, j) = EDGE_COLOR;
            else if(index > 0 && index <= seedCount)
                clredImage.at<Vec3b>(i, j) = colorTab[index - 1];
            else
                clredImage.at<Vec3b>(i, j) = DEFAULT_COLOR;
        }
    }
    delete[] colorTab;
}

int taskAMain(Mat & sMkrImage) {
    //initialize graphs
    Mat srcImage = imread("Resources/lena.png", 1);
    Mat uMkrImage(srcImage.size(), CV_8UC1, Scalar::all(0));
    Mat grayImage, wsImage(srcImage.size(), CV_8UC3);
    cvtColor(srcImage, sMkrImage, COLOR_BGR2GRAY);
    cvtColor(sMkrImage, grayImage, COLOR_GRAY2BGR);
    Mat temp;
    sMkrImage = Mat(srcImage.size(), CV_32S, Scalar::all(0));
    
    //input & calculate basic properties
    cout << "[WATERSHED] Enter number of seeds for watershed algorithm: ";
    cin >> seedCount;
    gridLength = sqrt(IMG_LENGTH * IMG_WIDTH / seedCount);
    minDist = gridLength;
    gridTtlCol = (int)((float)IMG_LENGTH / gridLength) + 1;
    gridTtlRow = (int)((float)IMG_WIDTH / gridLength) + 1;
    
    //generate ramdom seeds that fits requirements
    srand(int(time(0)));
    double dTime = (double)getTickCount();
    Seed * seedGroup = generateSeedGroup();
    dTime = (double)getTickCount() - dTime;
    cout << "[PHYSICS ENGINE] Seed generation completed in " << dTime / getTickFrequency() << " seconds using physics engine algorism." << endl;
    
    //watershed & show image
    markSeedOnImage(seedGroup, sMkrImage, MARK_AS_DOT);
    watershed(srcImage, sMkrImage);
    fillRandomColor(sMkrImage, wsImage);
    wsImage = grayImage * WATERSHED_IMG_TRANSPARENCY + wsImage * (1 - WATERSHED_IMG_TRANSPARENCY);
    markSeedOnImage(seedGroup, wsImage, MARK_AS_CIRCLE);
    imshow("taskA", wsImage);
    waitKey(0);
    delete[] seedGroup;
    return 0;
}
