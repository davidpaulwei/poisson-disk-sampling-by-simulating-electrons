#include "taskA.hpp"
#include "taskB.hpp"
#include "taskC.hpp"

int main() {
    Mat marker;
    taskAMain(marker);
    taskBMain(marker);
    taskCMain(marker);
    return 0;
}
