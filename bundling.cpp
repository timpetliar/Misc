//#include <stdio.h>      // fior printf()

using namespace std;
#include <iostream>     // for cpp std
#include <fstream>      // for ifstream, ofstream
#include <cmath>        // for ceil()

class DVD_Bundling
{
   #define GIDX(i,j) i*num_boxes+j      //index in array g
   #define NEW_RELEASE_WEEK 156
   #define WEEKS_TO_BUNDLING	8       //week before bundling occurs
   #define INFTY	1000
   #define BUNDLING_WEEK NEW_RELEASE_WEEK + WEEKS_TO_BUNDLING
//public data
public:

//private date
private:
    float eta;          // population size
    float phi;          // decision probability
    float nu;           // population "death" rate
    float psi;          // highest valuation in $
    int tau;            // week after which obsolescence begin
    float omega;        // limit of valuation discount rate
    int num_boxes;      
    float *prices;      // price data
    int prices_size;    // size of prices array

    float step_size;        
    float f;            // uniform density 
    float *g;           // The number of people with a valuations in intervals corresponding to coarsness of mesh
                        // it's a variable size array
    float *q;           // variable size array
    int q_size;         // size of q array                        
    float *revenue;
    int revenue_size;

    bool bundle;
    static const int prn_log = 1;

//public methods
public:
    DVD_Bundling(float _eta, float _phi, float _nu, float _psi, int _tau, float _omega, int _num_boxes, float *_prices, int _prices_size);
    ~DVD_Bundling();
    void simulate_path(int time);           //The main to generate model results
    float* get_rslt() { return q; }
    float* get_revenue(){return revenue;}
    /* print price array, e.g. for testing */
    void print_price_array()
    {
        cout << "Printing prices, size = " << prices_size << endl;
        int i;
        for (i=0; i<prices_size; i++) {
            cout << prices[i] << endl;
        }
    } 

//private methods
private:
    float __g_update__(float p_t_1, float p_t_2, float  a_t_1, float a_t_2);            //This method computes the new g(v)
    float __compute_q__(float integrate_from);          //This function approximated the integration used to obtain the quanity purchased
    float __discount__(int t);                          //This function computes the discount rate for each time value
    float __integrate_from__(float p_t, float a_t);     //lower bound for 'integration'
};

/* constructor */
DVD_Bundling::DVD_Bundling(float _eta, float _phi, float _nu, float _psi, int _tau, float _omega, int _num_boxes, float *_prices, int _prices_size)
{
   int i, j; 
   if (prn_log)
        cout << "DVD_Bundling::DVD_Bundling()" << endl;

    eta = _eta;
    phi = _phi;
    nu = _nu;
    psi = _psi;
    tau = _tau;
    omega = _omega;
    num_boxes = _num_boxes;
    prices = _prices;
    prices_size = _prices_size;
    bundle = false;
    step_size = psi/(num_boxes - 1);
    f = (float)1/(psi*psi);      //uniform density
	
    g = new float[num_boxes * num_boxes]; 

    for (i=0; i < num_boxes; i++) {
        for (j=0; j < num_boxes; j++) {
            g[GIDX(i,j)] = eta * f * step_size*step_size;         // number of people in each box of the approxmation for the g(v)
        }
    } 
    g[GIDX(num_boxes-1,num_boxes-1)] = 0;                      //right bottom element
 
    q = new float[_prices_size];
    q_size = _prices_size;
    revenue = new float[_prices_size];
    revenue_size = _prices_size; 
}

/* desctructor */
DVD_Bundling::~DVD_Bundling()
{
    if (prn_log)
        cout << "DVD_Bundling::~DVD_Bundling()" << endl;

    delete[] g;
    delete[] q;
    delete[] revenue;
}

/* The main to generate model results */
void DVD_Bundling::simulate_path(int time)
{
    int t;
    float p_t_1, p_t_2, a_t_1, a_t_2, q_t, integrate_from;

    for (t=0; t<NEW_RELEASE_WEEK; t++) {  // time zero exclued
        p_t_1 = prices[t];
	p_t_2 = INFTY;
        a_t_1 = __discount__(t);
	a_t_2 = 1;
       // integrate_from = __integrate_from__(p_t, a_t);
        revenue[t] = __g_update__(p_t_1, p_t_2, a_t_1, a_t_2);
        //q_t = __compute_q__(integrate_from);
       // q[t] = q_t;
    }
    for (t; t < BUNDLING_WEEK; t++) {
	p_t_1 = prices[t];
	p_t_2 = prices[t - NEW_RELEASE_WEEK];
 	a_t_1 = __discount__(t);
	a_t_2 = __discount__(t - NEW_RELEASE_WEEK);
	revenue[t] = __g_update__(p_t_1, p_t_2, a_t_1, a_t_2);
    }
    for (t; t < time; t++){
	bundle = true;
	p_t_1 = prices[t];
    	p_t_2 = prices[t - NEW_RELEASE_WEEK];
	a_t_1 = __discount__(t);
	a_t_2 = __discount__(t - NEW_RELEASE_WEEK);
	revenue[t] = __g_update__(p_t_1, p_t_2, a_t_1, a_t_2);
    }
    return ;
}

/* */
float  DVD_Bundling::__g_update__(float p_t_1, float p_t_2, float  a_t_1, float a_t_2)

{
    int i;
    int j;
    float g_new;
    float p_1 = a_t_1 * p_t_1;
    float p_2 = a_t_2 * p_t_2;
    float q_t = 0;
    float rev = 0;
    float p_b = p_t_1 + p_t_2 - 5;
    for (i = 0; i < num_boxes; i++){
	    for (j = 0; j < num_boxes; j++){
		    if ((p_b <= (a_t_1*i + a_t_2 * j)* step_size) && ( bundle == true) ){
			g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;
			g[GIDX(0,0)] += phi * g[GIDX(i,j)];
			g[GIDX(i,j)] = g_new;
			rev += p_b * phi * g[GIDX(i,j)];
		}
		else if (p_t_1 <= a_t_1 * i * step_size){
			g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;  
	   		g[GIDX(0,j)] += phi * g[GIDX(i,j)];
		   	g[GIDX(i,j)] = g_new;
		 rev += p_t_1 * phi * g[GIDX(i,j)];
		}
		else if(p_t_2 <= a_t_2  * j * step_size) {
			g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;
		        g[GIDX(i,0)] += phi * g[GIDX(i,j)];
			g[GIDX(i,j)] = g_new;
			rev += p_t_2 * phi * g[GIDX(i,j)]; 
		}
		else {
			g_new =(1- nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;
			g[GIDX(i,j)] = g_new;
		}

	    }
	}
    
    /*  else{
    for (i=0; i<num_boxes; i++)
    {
        if (i<integrate_from) {
            g_new =(1- nu) * g[i] + nu * eta * f * step_size * step_size;
            g[i] = g_new;
        }
        else {
            g_new = (1 - phi) *(1 - nu) * g[i] + nu * eta * f * step_size * step_size;
            g[0] += (1 - phi) *(1 - nu) * g[i]; //Those who have already purchased go to zero bucket
            g[i] = g_new;
        } 
    } */

    return rev;
}

/* */
float DVD_Bundling::__compute_q__(float integrate_from)
{
    int i;
    float q_t = 0.0;     //initialize to zero (no buyers) 
    
    if (integrate_from > num_boxes - 1) {
    }
    else {
        for (i=(int) integrate_from; i<num_boxes; i++) {
            q_t += g[i];    //sum of conusmer that have valuation greater than bound for summation (appromxation to integration)
        }
    }
    return (q_t * phi);     //multiply final resuly by phi, the awareness/probability of purchase parameter
}

/*  */
float DVD_Bundling::__discount__(int t)
{
    float a_t;
    if (t <= tau)
        a_t = 1;
    else
       a_t = omega + (1-omega) * ( (float)tau/(float)t ); 
    return a_t;
}

/* */
float DVD_Bundling::__integrate_from__(float p_t, float a_t) 
{
    float integrate_from = 0;
    // converts from indexing by value to indexing by box number 
    integrate_from = ceil(p_t/(a_t * step_size));
    return integrate_from;
}


#define TSIZE 260   //table size - num lines in the input file

int main()
{
    float t_val[TSIZE];
    float q_val[TSIZE];
    float p_val[TSIZE];

    ifstream ifs;
    ofstream ofs;

    int i=0, t=0;
    char c;
    string line; 

    cout << "Program starts" << endl;

    //open input & otput files
    ifs.open("pqdata.csv");
    if (! ifs.is_open() ) {
        cout << "Error. Cannot open input file" << endl;  
        return 4; 
    }

    ofs.open("bundling.csv");
    if (! ofs.is_open() ) {
        cout << "Error. Cannot open output file" << endl;
        return 4;
    }
 
    // read input data
    getline(ifs, line);         //skip 1st line
    //cout << "line " << line << endl;
    while ((ifs >> t >> c >> p_val[i] >> c >> q_val[i]) && (c==',')) { i++; }

    //example from python
    //test = DVD_Bundling(2000, .03,.005, 50, 3,  0.20, 501, p_val)
    //test.simulate_path(260)

    DVD_Bundling test(2000, .03, .005, 50, 3,  0.20, 501, p_val, TSIZE);
    //test.print_price_array(); 
    test.simulate_path(TSIZE);

    //put result data to file
    float *revenue = test.get_revenue();
    ofs << "t" << "," << "revenue" << endl;    
    for (t=0; t<TSIZE; t++) {
        ofs << t << "," << revenue[t] << endl;    
    }


    // close open files
    ifs.close();
    ofs.close();
    cout << "Program ends" << endl;
}
