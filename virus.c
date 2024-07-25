#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

int uniformeInt(float a, float b)
{   
    return a + (b - a) * genrand_real1();
}

float uniformFloat(float a, float b)
{   
    return a + (b - a) * genrand_real1();
}


#define GRID_SIZE 50
#define NUM_PEOPLE 100
#define NUM_DAYS 180 // 6 mois

#define HEALTHY 0
#define INCUBATING 1
#define INFECTIOUS 2
#define RECOVERED 3

typedef struct {
    int x;
    int y;
    int state;
    int day_infected;
} Person;

Person people[NUM_PEOPLE];
Person* grid[GRID_SIZE][GRID_SIZE];

void initialize_grid(int num_infected) {
    srand(time(NULL));

    // Initialize all grid cells to NULL
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            grid[i][j] = NULL;
        }
    }

    // Initialize people at random positions
    for (int i = 0; i < NUM_PEOPLE; i++) {
        int x = (int)(uniformFloat(0, 1) * GRID_SIZE);
        int y = (int)(uniformFloat(0, 1) * GRID_SIZE);

        // If the position is already occupied, find another one
        while (grid[x][y] != NULL) {
            x = (int)(uniformFloat(0, 1) * GRID_SIZE);
            y = (int)(uniformFloat(0, 1) * GRID_SIZE);
        }

        people[i].x = x;
        people[i].y = y;
        people[i].state = i < num_infected ? INCUBATING : HEALTHY;
        people[i].day_infected = i < num_infected ? 0 : -1;

        grid[x][y] = &people[i];
    }
}

void move_people() {
    for (int i = 0; i < NUM_PEOPLE; i++) {
        int old_x = people[i].x;
        int old_y = people[i].y;

        // Move the person in a random direction
        int new_x = (old_x + (rand() % 3 - 1) + GRID_SIZE) % GRID_SIZE;
        int new_y = (old_y + (rand() % 3 - 1) + GRID_SIZE) % GRID_SIZE;

        // If the new position is already occupied, don't move
        if (grid[new_x][new_y] != NULL) {
            continue;
        }

        people[i].x = new_x;
        people[i].y = new_y;

        // Update the grid
        grid[old_x][old_y] = NULL;
        grid[new_x][new_y] = &people[i];
    }
}

void update_states() {
    for (int i = 0; i < NUM_PEOPLE; i++) {
        if (people[i].state == INCUBATING && people[i].day_infected >= 2) {
            people[i].state = INFECTIOUS;
        } else if (people[i].state == INFECTIOUS && people[i].day_infected >= 12) {
            people[i].state = RECOVERED;
        }
        if (people[i].state == INCUBATING || people[i].state == INFECTIOUS) {
            people[i].day_infected++;
        }
    }
}

void collect_statistics(int* num_healthy, int* num_incubating, int* num_infectious, int* num_recovered) {
    *num_healthy = 0;
    *num_incubating = 0;
    *num_infectious = 0;
    *num_recovered = 0;
    for (int i = 0; i < NUM_PEOPLE; i++) {
        switch (people[i].state) {
            case HEALTHY:
                (*num_healthy)++;
                break;
            case INCUBATING:
                (*num_incubating)++;
                break;
            case INFECTIOUS:
                (*num_infectious)++;
                break;
            case RECOVERED:
                (*num_recovered)++;
                break;
        }
    }
}

float get_infection_probability(int day_infected) {
    if (day_infected == 2) {
        return 0.6;
    } else if (day_infected == 3) {
        return 0.8;
    } else if (day_infected == 4) {
        return 0.7;
    } else if (day_infected == 5) {
        return 0.6;
    } else if (day_infected == 6) {
        return 0.5;
    } else if (day_infected == 7) {
        return 0.4;
    } else if (day_infected == 8) {
        return 0.3;
    } else if (day_infected == 9) {
        return 0.2;
    } else if (day_infected == 10) {
        return 0.1;
    } else {
        return 0.0;
    }
}

void propagate_virus() {
    for (int i = 0; i < NUM_PEOPLE; i++) {
        if (people[i].state != INFECTIOUS) {
            continue;
        }

        // Check the 8 surrounding cells
        for (int dx = -1; dx <= 1; dx++) {
            for (int dy = -1; dy <= 1; dy++) {
                int x = (people[i].x + dx + GRID_SIZE) % GRID_SIZE;
                int y = (people[i].y + dy + GRID_SIZE) % GRID_SIZE;

                // If the cell is occupied by a healthy person, infect them with a certain probability
                if (grid[x][y] != NULL && grid[x][y]->state == HEALTHY && uniformFloat(0, 1) < get_infection_probability(people[i].day_infected)) {
                    grid[x][y]->state = INCUBATING;
                    grid[x][y]->day_infected = 0; // Start counting the days of infection
                }
            }
        }
    }
}

int report_results() {
    int num_healthy, num_incubating, num_infectious, num_recovered;
    collect_statistics(&num_healthy, &num_incubating, &num_infectious, &num_recovered);
    // printf("Healthy: %d\n", num_healthy);
    // printf("Incubating: %d\n", num_incubating);
    // printf("Infectious: %d\n", num_infectious);
    // printf("Recovered: %d\n", num_recovered);
    return num_infectious + num_recovered; // Return the total number of infected people
}

void calculate_confidence_interval(float* data, int size, float* mean, float* lower, float* upper) {
    // Calculate the mean
    *mean = 0;
    for (int i = 0; i < size; i++) {
        *mean += data[i];
    }
    *mean /= size;

    // Calculate the standard deviation
    float std_dev = 0;
    for (int i = 0; i < size; i++) {
        std_dev += pow(data[i] - *mean, 2);
    }
    std_dev = sqrt(std_dev / size);

    // Calculate the standard error
    float std_err = std_dev / sqrt(size);

    // Calculate the 95% confidence interval
    *lower = *mean - 1.96 * std_err;
    *upper = *mean + 1.96 * std_err;
}

int main() {
    int num_infected_start[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int num_days[] = {30, 90, 180, 365}; // 1 month, 3 months, 6 months, 1 year
    float results[sizeof(num_infected_start) / sizeof(int) * sizeof(num_days) / sizeof(int)];

    int k = 0;
    for (int i = 0; i < sizeof(num_infected_start) / sizeof(int); i++) {
        for (int j = 0; j < sizeof(num_days) / sizeof(int); j++) {
            initialize_grid(num_infected_start[i]);
            for (int day = 0; day < num_days[j]; day++) {
                move_people();
                propagate_virus();
                update_states();
            }
            results[k++] = report_results();
        }
    }

    float mean, lower, upper;
    calculate_confidence_interval(results, k, &mean, &lower, &upper);
    printf("Mean: %.2f\n", mean);
    printf("95%% Confidence interval: [%.2f, %.2f]\n", lower, upper);

    return 0;
}
