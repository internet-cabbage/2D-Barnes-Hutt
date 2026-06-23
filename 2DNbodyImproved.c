#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>




// ================================================================
// Time testing
// ================================================================

// variables used to measure time to execute certain sections
double buildTime = 0.0;
double forceTime = 0.0;
double accelTime = 0.0;

void printProgress(int stepVal, int tSteps, int width) {
    double perVal = (double) stepVal / (tSteps) * 100;

    char progBar[] = "========================================================================================================================";
    char emptyBar[] = "------------------------------------------------------------------------------------------------------------------------";
    char empty[] = "";

    // Width of the filled bar
    int filledWidth = (int) (perVal / 100 * width);
    int emptyPad = width - filledWidth;
    printf("%.s", emptyBar);
    // Bunch of mumbo jumbo I took ages to figure out
    // First bit (%5.2f%) just means 'print this number to 2 decimal places, and make it take up at least 5 places'
    // %.*s Basically just means 'Input how many characters to print of a string, input the string to be truncated'
    printf("\r %5.2f%% [%.*s%.*s]", perVal, filledWidth, progBar, emptyPad, emptyBar);
    fflush(stdout);
}

// ================================================================
// Structure definition
// ================================================================

typedef struct {
    double x,y;
} vec2;

typedef struct {
    vec2 position, velocity;
    double mass;
    // Stores a unique identifier for each body, to avoid self comparison
    int id;
} body;

typedef struct node {
    // Coordinates of the bounding box of the node
    double xmin, xmax, ymin, ymax;
    // Total mass of the node
    double massTot;
    // Position of the center of mass
    double cx,cy;



    // If the node is a leaf node, this variable stores the body data
    body nodeBody;
    // the children in the node [top left, top right, bottom left, bottom right]
    /*
    [0,1]
    [2,3] 
    */
    struct node* children[4];

} node;

// ================================================================
// Memory management function
// ================================================================

/*
Instead of storing all the node data in different places in cache, I can allocate one large pool of memory and store them all contiguously.
*/
node *pool; // A pointer to where all the nodes will be stored
int nodeCap; // The max amount of nodes in the pool
int nodeCount; // How many nodes are currently in the pool

// A bump allocator used to simplify the process of freeing all the memory
node* poolAlloc() {
    if (nodeCount >= nodeCap) {
        printf("Pool overflow error \n");
    }
    // The adress at which the new node is to be stored
    node* nodePtr = &pool[nodeCount];
    nodeCount++;
    *nodePtr = (node){0}; // Clears the node of previous data

    // Returns a pointer to where the node should be stored
    return nodePtr;
}

// ================================================================
// Quadtree helper functions
// ================================================================

// This function returns 1 if the body is within the node, and 0 otherwise
inline int isBodyInside(node *tNode, body *b) {
    // tNode just means targetNode, but it takes up less space on my screen
    if ((tNode->xmin < b->position.x) && (b->position.x < tNode->xmax) && (tNode->ymin < b->position.y) && (b->position.y < tNode->ymax)) {
        //printf("Node x-range (%f <-> %f), y-range (%f <-> %f) \nBody Coords (%f, %f) \nBody id: (ID %d) \n \n", tNode->xmin, tNode->xmax, tNode->ymin, tNode->ymax, b->position.x, b->position.y, b->id);
        return 1;
    }
    else {
        return 0;
    }
}

// This returns the length of the sides of the nodes, assuming that the node is a square

inline double lenth(node *tNode) {
    return (tNode->xmax - tNode->xmin);
}

void quadDivide(node *tNode) {
    double xmid = (tNode->xmin + tNode->xmax) / 2;
    double ymid = (tNode->ymin + tNode->ymax) / 2;

    /*
    [0,1]
    [2,3] 
    */

    node emptyNode;

    tNode->children[0] = poolAlloc();
    *tNode->children[0] = (node){.xmin = tNode->xmin, .xmax = xmid, .ymin = ymid, .ymax = tNode->ymax};

    tNode->children[1] = poolAlloc();
    *tNode->children[1] = (node){.xmin = xmid, .xmax = tNode->xmax, .ymin = ymid, .ymax = tNode->ymax};
    
    tNode->children[2] = poolAlloc();
    *tNode->children[2] = (node){.xmin = tNode->xmin, .xmax = xmid, .ymin = tNode->ymin, .ymax = ymid};
    
    tNode->children[3] = poolAlloc();
    *tNode->children[3] = (node){.xmin = xmid, .xmax = tNode->xmax, .ymin = tNode->ymin, .ymax = ymid};

}

// This little function serves to help the insertBody function, by calculating
// which index position to insert the new node at
inline int insertIndex(node *parentNode, double x, double y) {
    double xmid = (parentNode->xmin + parentNode->xmax) / 2.0;
    double ymid = (parentNode->ymin + parentNode->ymax) / 2.0;
    
    // Top nodes
    if (ymid < y) {
        // Is node on the right?
        if (xmid < x) {
            return 1;
        }
        // top left node
        else {
            return 0;
        }
    }
    // bottom nodes
    else {
        // bottom right node
        if (xmid < x) {
            return 3;
        }
        // bottom left node
        else {
            return 2;
        }
    }

}

void insertBody(node *targetNode, body *b) {

    /* The process for inserting goes as follows:
        1) Check if the body being inserted, is even within the node
        2) if so, check if the node has children
            3.a) if yes, then recursively insert this body into a child node
            3.b) if no, then check if this node is occupied:
                4.a) If the node is occupied (i.e a leaf node), give this node 4 kids, and insert the old node
                into one of the child nodes, and do the same for the node being inserted. If they both get inserted into the same node
                then recursively repeat this process until they are at different locations
                4.b) If the node isnt occupied, set the node info to that of the body
    */

    // If it returns one, then the body IS in the node

    if (isBodyInside(targetNode, b) == 1) {
        // returns one if the targetNode does have kids
        // If it has kids, then we update this node, and insert the body into one of its kids
        if (targetNode->children[0] != NULL) {

            // Firstly we have to update the center of mass and total mass of the parent node

            double newTotalMass = targetNode->massTot + b->mass;
            
            // This is just the standard centre of mass formula from physics
            targetNode->cx = ((targetNode->cx * targetNode->massTot)+(b->position.x * b->mass)) / newTotalMass;
            targetNode->cy = ((targetNode->cy * targetNode->massTot)+(b->position.y * b->mass)) / newTotalMass;
            targetNode->massTot = newTotalMass;
            
            // Finds the index position of the node the body should be inserted in
            int index = insertIndex(targetNode, b->position.x, b->position.y);
            // Recursively calls the function until it is inserted
            insertBody(targetNode->children[index],b);
        }
        else {
            // Check if the node is a leaf node
            if (targetNode->massTot != 0.0) {
                // Finds the location to put the 'native' node and new node into
                int nativeIndex = insertIndex(targetNode, targetNode->nodeBody.position.x, targetNode->nodeBody.position.y);
                int newIndex = insertIndex(targetNode, b->position.x, b->position.y);

                // Creates 4 child nodes
                quadDivide(targetNode);

                // The native node gets inserted into the correct child, then removed from the native node
                insertBody(targetNode->children[nativeIndex],&targetNode->nodeBody);
                targetNode->nodeBody = (body) {0.0};
                targetNode->nodeBody.id = -1;

                // Inserts the new node into the correct place
                insertBody(targetNode->children[newIndex], b);

                // Now we have to recompute the centre of mass and its location
                double newTotalMass = targetNode->massTot + b->mass;
            
                targetNode->cx = ((targetNode->cx * targetNode->massTot)+(b->position.x * b->mass)) / newTotalMass;
                targetNode->cy = ((targetNode->cy * targetNode->massTot)+(b->position.y * b->mass)) / newTotalMass;
                targetNode->massTot = newTotalMass;
                
            }
            else {
                // If the targetNode does not have kids, then it is a leaf node, so we insert the body here
                // This updates the value of the 'nodeBody' attribute of the targetNode struct, to be equal to the value held by pointer b
                targetNode->nodeBody = *b;
                targetNode->massTot = b->mass;
                targetNode->cx = b->position.x;
                targetNode->cy = b->position.y;
            }
        }
    }

    /*
    // If it returns zero, then the body wasnt inside the node
    else if (isBodyInside(targetNode, b) == 0) {
        printf("Body (ID: %d) being inserted, isn't in node being inserted to. \n", b->id);
    }
    */
    
    // If the node is empty then nodeBody == NULL
}
// ================================================================
// Force calculator function
// ================================================================

// This calculates the force that a specific node 'tNode' is acting on the body 'b'
vec2 calculateForce(node *tNode, body *b, double antiSingularity, double G, double Theta) {

    // If its total mass is zero, it contains no bodies and thus exerts no force
    if (tNode->massTot == 0.0) {
        return (vec2) {0.0,0.0};
    }
    /* if the node being considered is the one containing the body, then it doesnt exactly exert any force on itself does it?
     The body is the same as the one within the node if they have the same pointer 
     you cant directly compare structs in c, so you gotta do this for that reason as well

     THE ABOVE IS A LIE!!! I am leaving it there incase I make the same error. But the node only stores a copy of the object, not its pointer
     so the body in the node and the real body will have different adresses
    */
    else if (tNode->nodeBody.id == b->id) {
        return (vec2) {0.0,0.0};
    }

    // Calculates the vector distance between the center of mass of the node, and the body position
    else {
        double dx = tNode->cx - b->position.x;
        double dy = tNode->cy - b->position.y;
        double r2 = (dx * dx) + (dy * dy) + (antiSingularity * antiSingularity);
        

        double nodeSize = lenth(tNode);

        // This decides whether or not to use the approximation.
        // If the node has no child nodes, applying the approximation is the same as directly summing it
        // Instead of the standard nodeSize/r < theta comparison, I squared it so that I dont have to perform an unnecesary square root
        if (tNode->children[0] == NULL || ((nodeSize*nodeSize)/r2) < Theta*Theta) {
            
            double fMag = (G * tNode->massTot * b->mass) / (r2);
            double r = sqrt(r2);
            double fx = fMag * (dx/r);
            double fy = fMag * (dy/r);
            return (vec2) {fx,fy};
        }
        // If the nodes are too close to approximate, it recursively calculates it directly
        else {
            vec2 totalForce = {0.0,0.0};
            for (int i = 0; i < 4; i++) {
                if (tNode->children[i] != NULL) {
                    vec2 recursiveForce = calculateForce(tNode->children[i], b, antiSingularity, G, Theta);
                    totalForce.x += recursiveForce.x;
                    totalForce.y += recursiveForce.y;
                }
            }
            return totalForce;
        }
    }
}





// ================================================================
// Quadtree helper functions
// ================================================================

double* randomGen(int lower, int upper, unsigned int N) {
    // Adress at which the array is saved at
    double *adress = calloc(N, sizeof(double));

    for (int i = 0; i < N; i++) {
        double val = (rand() % (upper-lower+1)) + lower;
        *(adress + i) = val;
    }
    return adress;
}

double* randomContinuous(double maxMag, int N) {
    double *adress = calloc(N, sizeof(double));

    for (int i = 0; i < N; i++) {
        if ((double)rand()/RAND_MAX > 0.5) {
            double val = maxMag * rand()/RAND_MAX;
            *(adress + i) = val;
        }
        else {
            double val = - maxMag * rand()/RAND_MAX;
            *(adress + i) = val;
        }
    }
    return adress;
}

// Cute little helper function to make it easier to write data to the binary file
void writeFrame(FILE *dataFile, body *bodies, unsigned int N, float *frameBuffer) {
    for (int i = 0; i < N; i++) {
        frameBuffer[2*i] = (float)bodies[i].position.x;
        frameBuffer[(2*i)+1] = (float)bodies[i].position.y;
    }
    fwrite(frameBuffer, sizeof(float), 2*N, dataFile);
}

void timeLoop(body *bodies, int N, double antiSingularity, double G, int tSteps, double dt, double theta, int xmax, int ymax) {
    // The array storing the acceleration info for the bodies
    vec2 *accelArray = calloc(N, sizeof(vec2));

    /*
    When outputed to the binary file, I need to specify the value of N,
    as it determines how large each frame is in the file. So I just output the value N
    as the first item in the file.
    */

    // So firstly I have to create the entire tree datastructure, by recursively inserting each node

    nodeCount = 0;
    node* rootPtr = poolAlloc();
    *rootPtr = (node){0};
    rootPtr->xmax = xmax * 2; rootPtr->xmin = -xmax * 2;
    rootPtr->ymax = ymax * 2; rootPtr->ymin = -ymax * 2;

    double t0 = omp_get_wtime();
    // Inserts all bodies into the tree to populate it
    for (int i = 0; i < N; i++) {
        insertBody(rootPtr, &bodies[i]);
    }
    
    double t1 = omp_get_wtime();
    // Uses the lovely tree to calculate the forces and accelerations on all the bodies
    // It also makes use of openmp to paralellise it
    #pragma omp parallel for schedule(dynamic, 32)
    for (int i = 0; i < N; i++) {
        vec2 f = calculateForce(rootPtr, &bodies[i], antiSingularity, G, theta);
        
        accelArray[i].x = f.x / bodies[i].mass;
        bodies[i].velocity.x += accelArray[i].x * dt;
        bodies[i].position.x += bodies[i].velocity.x * dt;


        accelArray[i].y = f.y / bodies[i].mass;
        bodies[i].velocity.y += accelArray[i].y * dt;
        bodies[i].position.y += bodies[i].velocity.y * dt;
    }
    double t2 = omp_get_wtime();


    buildTime += t1 - t0;
    forceTime += t2 - t1;

    free(accelArray);
    
    
}


int main() {
    fprintf(stderr,"\nMax threads: %d \n", omp_get_max_threads());
    // Seeding the random generation function, so it doesnt repeat values
    srand(time(NULL));
    // Simulation parameters
    unsigned int N = 20000; // Number of bodies
    int tSteps = 20000;
    double antiSingularity = 1.0;
    double G = 1;
    double theta = 1.0; 
    int xMax = 4000; // Maximum x distance
    int yMax = 4000; // Self explanatory
    int vMax = 66; // Max velocity
    int mMax = 100; // Max mass
    int mMin = 90;
    double dt = 0.04;
    
    double *xVals = randomContinuous(xMax, N);
    double *yVals = randomContinuous(yMax, N);

    double *vxVals = randomContinuous(vMax, N);
    double *vyVals = randomContinuous(vMax, N);
    double *mVals = randomGen(mMin,mMax, N);

    // Initialise bodies loop

    body *bodies = calloc(N, sizeof(body));
    for (int i = 0; i < N; i++) {
        vec2 position, velocity;
        position.x = xVals[i];
        position.y = yVals[i];

        velocity.x = vxVals[i];
        velocity.y = vyVals[i];

        bodies[i].position = position;
        bodies[i].velocity = velocity;
        bodies[i].mass = mVals[i]; 
        bodies[i].id = i;
    }

    free(xVals);
    free(yVals);
    free(vxVals);
    free(vyVals);
    free(mVals);
    xVals = yVals = vxVals = vyVals = mVals = NULL;

    /*
    In order for my chosen integrator to work, I need to calculate the initial acceleration, so I can 'offset' the velocity
    by half of a time step. This integrator method is known as the leapfrog integrator, and it is better at conserving energy than
    the euler step method I was originally going to use.    
    */ 

    // Using the pool memory allocator, for cache locality

    // nodeCap determines how many nodes can be stored in the memory pool, any more nodes than this will cause the pool to overflow and the program to presumably become very corrupted

    nodeCap = 8 * N;
    pool = calloc(nodeCap, sizeof(node));
    nodeCount = 0;

    if (pool == NULL) {
        printf("Pool memory allocation failed!!! PANIC!!! \n");
        exit(1);
    }

    fprintf(stderr, "Initial tree build successful. \n");
    // Creating the root node
    node* rootPtr = poolAlloc();
    *rootPtr = (node){0};
    rootPtr->xmax = xMax * 2;
    rootPtr->xmin = -xMax * 2;
    rootPtr->ymax = yMax * 2;
    rootPtr->ymin = -yMax * 2;
    // Creating the tree

    // Inserts all bodies into the tree to populate it
    for (int i = 0; i < N; i++) {
        insertBody(rootPtr, &bodies[i]);
    }
    // Uses the lovely tree to calculate the forces and accelerations on all the bodies
    vec2 *accelArray = calloc(N, sizeof(vec2));
    for (int i = 0; i < N; i++) {
        vec2 f = calculateForce(rootPtr, &bodies[i], antiSingularity, G, theta);
        accelArray[i].x = f.x / bodies[i].mass;
        accelArray[i].y = f.y / bodies[i].mass;
    }

    // Offsets the velocity to apply the leapfrog correction

    for (int i = 0; i < N; i++) {
        bodies[i].velocity.x -= accelArray[i].x * (dt/2); 
        bodies[i].velocity.y -= accelArray[i].y * (dt/2); 
    }
    
    // Now I free all the memory used to perform those calculations

    //freeTree(rootPtr);
    free(accelArray);

    fprintf(stderr, "Initial tree collapse successful. \n");

    // This is the main loop which performs all the calculations, and saving the data

    // File to write to
    FILE *dataFile = fopen("DataOutput.bin", "wb");

    // Checking if the file was actually created, as who knows what the C program would do otherwise
    if (dataFile == NULL) {
        fprintf(stderr ,"Output file can't be created:");
    }
    else {
        fprintf(stderr, "Output file successfully created. \n \nBeginning calculations: \n \n");
    }

    // Writes the parameters required for the PythonRenderer to interpret the data
    fwrite(&N,sizeof(N),1,dataFile);
    fwrite(&tSteps,sizeof(tSteps),1,dataFile);

    // The array storing the positions of all the bodies, to be written to a buffer
    float * frameBuffer = calloc(2*N, sizeof(float));

    // Initial conditions of the system
    writeFrame(dataFile,bodies,N,frameBuffer);

    // Code for the loading bar
    static char barString[] = "||||||||||||||||||||||||||||||||||||||||";
    static int barLength = 40;

    int barInterval = (tSteps / (6 * barLength)) ;

    for (int j = 0; j < tSteps; j++) {
        if ((j+1) % barInterval == 0) {
            printProgress(j, tSteps, barLength);
        }
        timeLoop(bodies,N,antiSingularity,G,tSteps,dt,theta,xMax,yMax);
         /* And finally writing this output to a file, each row will contain the state for a single body at a single time
        [xPos_0, yPos_0, xPos_1, yPos_1...]
        */
        writeFrame(dataFile,bodies,N,frameBuffer);
    }
    fclose(dataFile);
    free(frameBuffer);
    free(bodies);
    free(pool);

    fprintf(stderr, "\nBuild time: %.2f, force time: %.2f \n", buildTime, forceTime);
}