using namespace std;
#include <iostream>     // for cpp std
#include <fstream>      // for ifstream, ofstream
#include <cmath>        // for ceil()

class DVD_Bundling
{
   #define GIDX(i,j) i*num_boxes+j      //index in array g
   #define NEW_RELEASE_WEEK 156
   #define WEEKS_TO_BUNDLING	52       //week before bundling occurs
   #define INFTY	1000
   #define BUNDLING_WEEK NEW_RELEASE_WEEK + WEEKS_TO_BUNDLING	//week at which bundling begins
   #define TSIZE 260   //table size - num lines in the input file

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
    int bundling_week;  // /week at which bundling begins

    float step_size;        
    float f;            // uniform density 
    float *g;           // The number of people with a valuations in intervals corresponding to coarsness of mesh
                        // it's a variable size array
                            
    float *revenue;	// pointer to array storing weekly revenue
    int revenue_size;	//number of weeks in revenue array
    
    bool bBundle;	// true if bundling, false if no bundling
    static const int prn_log = 0;

//public methods
public:
    DVD_Bundling(float _eta, float _phi, float _nu, float _psi, int _tau, float _omega, 
                 int _num_boxes, float *_prices, int _prices_size, int _bundling_week);
    ~DVD_Bundling();
    void simulate_path(int time);           //The main to generate model results
    float* get_revenue(){return revenue;}   // allows use of revenue array outside class
    float  comp_NPV();	 		    //Compute Net Present Value at t = 0	
    float buy_bundle(float p, int i, int j);
    float buy_1(float p, int i, int j);
    float buy_2(float p, int i, int j);
    void  buy_none(int i, int j);   
    float buy_both(float p1, float p2, int i, int j);

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
    float g_update(float p_t_1, float p_t_2, float  a_t_1, float a_t_2);	//This method computes the new g(v) and then computes weekly sales
    float discount(int t);                         				//This function computes the discount rate for each time value
};

/* constructor */
DVD_Bundling::DVD_Bundling(float _eta, float _phi, float _nu, float _psi, int _tau, float _omega, 
                           int _num_boxes, float *_prices, int _prices_size, int _bundling_week)
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
    bundling_week = _bundling_week;
    bBundle = false;
    step_size = psi/(num_boxes - 1);
    f = (float)1/(psi*psi);      //uniform density
	
    g = new float[num_boxes * num_boxes]; 

    for (i=0; i < num_boxes; i++) {
        for (j=0; j < num_boxes; j++) {
			if( (i == num_boxes -1) || (j == num_boxes -1) ) {
				g[GIDX(i,j)] = 0;
			}
			else{	
				g[GIDX(i,j)] = eta * f * step_size*step_size;
			}			// number of people in each box of the approxmation for the g(v)
        }
    } 
   //[GIDX(num_boxes-1,num_boxes-1)] = 0;                      //right bottom element
    revenue = new float[_prices_size]; 				//allocate memory for revenue array and assign to revenue pointer
    revenue_size = _prices_size; 				//size of revenue array
}

/* desctructor */
DVD_Bundling::~DVD_Bundling()
{
    if (prn_log)						//used for debugging
        cout << "DVD_Bundling::~DVD_Bundling()" << endl;

    delete[] g;							//free memory allocated for g array
    delete[] revenue;						//free memory allocated for revenue array
}

/* The main to generate model results */
void DVD_Bundling::simulate_path(int time)
{
    int t;
    float p_t_1, p_t_2, a_t_1, a_t_2; 			//weekly price and dicount factor for DVD1 and DVD2

    for (t=0; t<NEW_RELEASE_WEEK; t++) {			//Only first DVD has been released
        p_t_1 = prices[t];	
	p_t_2 = INFTY;						//this allows writing all function for 2 DVD case
        a_t_1 = discount(t);	
	a_t_2 = 1;						//arbitrary placeholder value
        revenue[t] = g_update(p_t_1, p_t_2, a_t_1, a_t_2); 	//update g and store the revenue for the week
    }
    for (; t < bundling_week; t++) {				//new DVD released. NO BUNDLING YET
	p_t_1 = prices[t];
	//p_t_2 = INFTY;
	p_t_2 = prices[t - NEW_RELEASE_WEEK];	
 	a_t_1 = discount(t);
	a_t_2 = discount(t - NEW_RELEASE_WEEK);
	revenue[t] = g_update(p_t_1, p_t_2, a_t_1, a_t_2);
    }
    for (; t < time; t++){					//Bundle offered. Continue until end of tracking period
	bBundle = true;		//bundle se calculation will now be activated
	p_t_1 = prices[t];
	//p_t_2 = INFTY; 
    	p_t_2 = prices[t - NEW_RELEASE_WEEK];
	a_t_1 = discount(t);
	a_t_2 = discount(t - NEW_RELEASE_WEEK);
	revenue[t] = g_update(p_t_1, p_t_2, a_t_1, a_t_2); 
    }
    return ;
}

/* */
float DVD_Bundling::buy_bundle(float p_b, int i, int j){
	float g_new;
	g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;
	g[GIDX(0,0)] += phi * g[GIDX(i,j)];							//Already bought both
	g[GIDX(i,j)] = g_new;
	return (p_b * phi * g[GIDX(i,j)]);
}

/* */
float DVD_Bundling::buy_1 (float p_t_1, int i,int j){
	float g_new;
	g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;  
	g[GIDX(0,j)] += phi * (1 - nu) * g[GIDX(i,j)];							//May buy DVD2  in future
	g[GIDX(i,j)] = g_new;
	return (p_t_1 * phi * g[GIDX(i,j)]);
}

float DVD_Bundling::buy_2 (float p_t_2,int i,int j){
	float g_new;
	g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;	
	g[GIDX(i,0)] += phi * (1 - nu) * g[GIDX(i,j)];							//May buy DVD1 in future
	g[GIDX(i,j)] = g_new;
	return  (p_t_2 * phi * g[GIDX(i,j)]);
}

float DVD_Bundling::buy_both(float p_t_1, float p_t_2, int i, int j){
	float g_new;
	g_new = (1- phi) * (1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size *step_size;
	g[GIDX(0,0)] += phi * (1 - nu) * g[GIDX(i,j)]; 
	g[GIDX(i,j)] = g_new;
	return ( (p_t_1 + p_t_2) * phi * g[GIDX(i,j)] );
}
void DVD_Bundling::buy_none (int i, int j){
	float g_new;
	g_new =(1- nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;
	g[GIDX(i,j)] = g_new;
	return;
} 
			
/* */
float  DVD_Bundling::g_update(float p_t_1, float p_t_2, float  a_t_1, float a_t_2)
{
    int i;
    int j;
    float rev = 0;
    //float p_b = INFTY;
    float p_b = p_t_1 + p_t_2 - 5;										//price of bundle
    for (i = 0; i < num_boxes; i++){
	    for (j = 0; j < num_boxes; j++){
		    if ( (p_b <= (a_t_1*(float)i + a_t_2 *(float) j)* step_size) && ( bBundle == true) ){			//Buy bundle ?
			
			if ( (a_t_1 * ((float) i * step_size)) < (p_t_1 - 5) ){
 				rev += buy_2(p_t_2, i, j);
			}
			if ( (a_t_2 *((float) j * step_size)) < (p_t_2 -5) ){
				rev += buy_1(p_t_1, i, j);
			}
			else {
				rev += buy_bundle(p_b, i, j);
			}

		   
			
		}
		else if (p_t_1 <= a_t_1 *(float) i * step_size){							//Buy only DVD1 ?
			if (p_t_2 <= a_t_2 *(float) j * step_size){							//Buy only DVD1 ?
				rev += buy_both(p_t_1, p_t_2, i, j);
			}
			else{
				rev += buy_1(p_t_1, i,j );
			}


		}
		else if(p_t_2 <= a_t_2  * (float) j * step_size) {
			rev += buy_2(p_t_2, i, j);									//But only DVD2 ?
			 
		}
		else {												//Buy nothing.
			buy_none(i, j);
			
		} 	

	    }			//close j loop
	}			//close i loop
    
   return rev;
}


/*  */
float DVD_Bundling::discount(int t)
{
    float a_t;
    if (t <= tau)
        a_t = 1;
    else
       a_t = omega + (1-omega) * ( (float)tau/(float)t ); 
    return a_t;
}


/* */
float DVD_Bundling::comp_NPV()
{
    int t; 
    float r = .02; //discount rate
    float npv = 0; //Net Present Value
    for (t = 0; t<TSIZE; t++){
	    npv += revenue[t] * (1/pow((1+r),(float) t));
    }
    //cout <<"\nNPV: " << npv << endl;
    return npv;
}


/* */
//#define TSIZE 260   //table size - num lines in the input file
int main()
{
    float q_val[TSIZE];
    float p_val[TSIZE];

    ifstream ifs;
    ofstream ofs;

    int i=0, t=0;
    char c;
    string line; 
    float npv;

    cout << "Program starts" << endl;

    //open input & otput files
    ifs.open("pqdata.csv");
    if (! ifs.is_open() ) {
        cout << "Error. Cannot open input file" << endl;  
        return 4; 
    }

    ofs.open("npv.csv");
    if (! ofs.is_open() ) {
        cout << "Error. Cannot open output file" << endl;
        return 4;
    }
 
    // read input data
    getline(ifs, line);         //skip 1st line
    //cout << "line " << line << endl;
    while ((ifs >> t >> c >> p_val[i] >> c >> q_val[i]) && (c==',')) { i++; }


    DVD_Bundling *pTest;    //pTest is pointer to instance of DVD_Bundling class
    
/*
    pTest = new DVD_Bundling(2000, .03, .005, 50, 3,  0.20, 501, p_val, TSIZE, BUNDLING_WEEK);
    //pTest->print_price_array(); 
    pTest->simulate_path(TSIZE);	//run the model
    npv = pTest->comp_NPV();
    cout <<"\nNPV: " << npv << endl;
    delete pTest;    
*/
     ofs << "t" << "," << "npv" << endl;    
 
    for (t=NEW_RELEASE_WEEK; t<TSIZE; t++) {
        pTest = new DVD_Bundling(2000, .03, .005, 50, 3,  0.20, 501, p_val, TSIZE, t);
        pTest->simulate_path(TSIZE);        //run the model
        npv = pTest->comp_NPV();
        cout << "t=" << t-NEW_RELEASE_WEEK << " NPV: " << npv << endl;    
        ofs << t-NEW_RELEASE_WEEK << "," << npv << endl;
        delete pTest; 
    }

/*
    //output results in csv file
    float *revenue = pTest->get_revenue();				
    ofs << "t" << "," << "revenue" << endl;    
    for (t=0; t<TSIZE; t++) {
        ofs << t << "," << revenue[t] << endl;    
    }
*/

    // close open files
    ifs.close();
    ofs.close();
    cout << "\nProgram ends" << endl;
}
