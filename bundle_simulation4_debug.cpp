using namespace std;
#include <iostream>     // for cpp std
#include <fstream>      // for ifstream, ofstream
#include <cmath>        // for ceil()

class DVD_Bundling
{
   #define GIDX(i,j) i*num_boxes+j      //index in array g
   #define NEW_RELEASE_WEEK 156
   #define WEEKS_TO_BUNDLING	52       //week before bundling occurs
   #define INFTY	10000
   #define BUNDLING_WEEK NEW_RELEASE_WEEK + WEEKS_TO_BUNDLING	//week at which bundling begins
   #define TSIZE 260   //table size - num lines in the input file

//public data
public:

//private date
private:
    double eta;          // population size
    double phi;          // decision probability
    double nu;           // population "death" rate
    double psi;          // highest valuation in $
    int tau;            // week after which obsolescence begin
    double omega;        // limit of valuation discount rate
    int num_boxes;      
    double *prices;      // price data
    int prices_size;    // size of prices array

    double step_size;        
    double f;            // uniform density 
    double *g;           // The number of people with a valuations in intervals corresponding to coarsness of mesh
                        // it's a variable size array
    
                            
    double *revenue;	// pointer to array storing weekly revenue
    int revenue_size;	//number of weeks in revenue array

    bool bundle;	// true if bundling, false if no bundling
    static const int prn_log = 1;

//public methods
public:
    DVD_Bundling(double _eta, double _phi, double _nu, double _psi, int _tau, double _omega, int _num_boxes, double *_prices, int _prices_size);
    ~DVD_Bundling();
    void simulate_path(int time);           //The main to generate model results
    double* get_revenue(){return revenue;}   // allows use of revenue array outside class
    void  NPV();	 		//Compute Net Present Value at t = 0	
    double buy_bundle(double p, int i, int j);
    double buy_1(double p, int i, int j);
    double buy_2(double p, int i, int j);
    void  buy_none(int i, int j);   
    double buy_both(double p1, double p2, int i, int j);

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
    double g_update(double p_t_1, double p_t_2, double  a_t_1, double a_t_2);	//This method computes the new g(v) and then computes weekly sales
    double discount(int t);                         				//This function computes the discount rate for each time value
};
/* constructor */
DVD_Bundling::DVD_Bundling(double _eta, double _phi, double _nu, double _psi, int _tau, double _omega, int _num_boxes, double *_prices, int _prices_size)
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
    f = (double)1/(psi*psi);      //uniform density
	
    g = new double[num_boxes * num_boxes]; 

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
    revenue = new double[_prices_size]; 				//allocate memory for revenue array and assign to revenue pointer
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
    double p_t_1, p_t_2, a_t_1, a_t_2; 			//weekly price and dicount factor for DVD1 and DVD2

    for (t=0; t<NEW_RELEASE_WEEK; t++) {			//Only first DVD has been released
        p_t_1 = prices[t];	
	p_t_2 = INFTY;						//this allows writing all function for 2 DVD case
        a_t_1 = discount(t);	
	a_t_2 = 1;						//arbitrary placeholder value
        revenue[t] = g_update(p_t_1, p_t_2, a_t_1, a_t_2); 	//update g and store the revenue for the week
    }
    for (; t < BUNDLING_WEEK; t++) {				//new DVD released. NO BUNDLING YET
	p_t_1 = prices[t];
	//p_t_2 = INFTY;
	p_t_2 = prices[t - NEW_RELEASE_WEEK];	
 	a_t_1 = discount(t);
	a_t_2 = discount(t - NEW_RELEASE_WEEK);
	revenue[t] = g_update(p_t_1, p_t_2, a_t_1, a_t_2);
    }
    for (; t < time; t++){					//Bundle offered. Continue until end of tracking period
	bundle = true;		//bundle se calculation will now be activated
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
double DVD_Bundling::buy_bundle(double p_b, int i, int j){
	double g_new;
	g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;
	g[GIDX(0,0)] += phi * g[GIDX(i,j)];							//Already bought both
	g[GIDX(i,j)] = g_new;
	return (p_b * phi * g[GIDX(i,j)]);
}

/* */
double DVD_Bundling::buy_1 (double p_t_1, int i,int j){
	double g_new;
	g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;  
	g[GIDX(0,j)] += phi * (1 - nu) * g[GIDX(i,j)];							//May buy DVD2  in future
	g[GIDX(i,j)] = g_new;
	return (p_t_1 * phi * g[GIDX(i,j)]);
}

double DVD_Bundling::buy_2 (double p_t_2,int i,int j){
	double g_new;
	g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;	
	g[GIDX(i,0)] += phi * (1 - nu) * g[GIDX(i,j)];							//May buy DVD1 in future
	g[GIDX(i,j)] = g_new;
	return  (p_t_2 * phi * g[GIDX(i,j)]);
}

double DVD_Bundling::buy_both(double p_t_1, double p_t_2, int i, int j){
	double g_new;
	g_new = (1- phi) * (1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size *step_size;
	g[GIDX(0,0)] += phi * (1 - nu) * g[GIDX(i,j)]; 
	g[GIDX(i,j)] = g_new;
	return ( (p_t_1 + p_t_2) * phi * g[GIDX(i,j)] );
}
void DVD_Bundling::buy_none (int i, int j){
	double g_new;
	g_new =(1- nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size;
	g[GIDX(i,j)] = g_new;
	return;
} 
			
/* */
double  DVD_Bundling::g_update(double p_t_1, double p_t_2, double  a_t_1, double a_t_2) {
    int i;
    int j;
    double rev = 0;
    //double p_b = INFTY;
    double p_b = p_t_1 + p_t_2 - 5;										//price of bundle
    for (i = 0; i < num_boxes; i++){
	    for (j = 0; j < num_boxes; j++){
		    if ( (p_b <= (a_t_1*(double)i + a_t_2 *(double) j)* step_size) && ( bundle == true) ){			//Buy bundle ?
			
			if ( (a_t_1 * (double) i * step_size) < (p_t_1 - 5) ){
 				rev += buy_2(p_t_2, i, j);
			}
			if( (a_t_2 *(double) j * step_size) < (p_t_2 -5) ){
				rev += buy_1(p_t_1, i, j);
			}
		/*	if (  ( (a_t_1 * (double) i * step_size) < (p_t_1 - 5) ) && ( (a_t_2 *(double) j * step_size) < (p_t_2 -5) ) ){
			;
			}*/ 
			//if (i == -5){;}
			//if (i ==  -10) {;}
			//if ( i == -6 ) {;}
			//if ( j == -5.0 ) {;}	
			//if ( 5 == 5) {}
			else {
				rev += buy_bundle(p_b, i, j);
			}

		   
			
		}
		else if (p_t_1 <= a_t_1 *(double) i * step_size){							//Buy only DVD1 ?
			if (p_t_2 <= a_t_2 *(double) j * step_size){							//Buy only DVD1 ?
				rev += buy_both(p_t_1, p_t_2, i, j);
			}
			if ( 5== -5) {;}
			else{
				rev += buy_1(p_t_1, i,j );
			}


		}
		else if(p_t_2 <= a_t_2  * (double) j * step_size) {
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
double DVD_Bundling::discount(int t)
{
    double a_t;
    if (t <= tau)
        a_t = 1;
    else
       a_t = omega + (1-omega) * ( (double)tau/(double)t ); 
    return a_t;
}


/* */
void DVD_Bundling::NPV()
{
    int t; 
    double r = .02; //discount rate
    double npv = 0; //Net Present Value
    for (t = 0; t<TSIZE; t++){
	    npv += revenue[t] * (1/pow((1+r),(double) t));
    }
    	cout <<"\nNPV: " << npv << endl;
	return;
}


/* */
//#define TSIZE 260   //table size - num lines in the input file
int main()
{
    double q_val[TSIZE];
    double p_val[TSIZE];

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

    ofs.open("bundling_revenue2.csv");
    if (! ofs.is_open() ) {
        cout << "Error. Cannot open output file" << endl;
        return 4;
    }
 
    // read input data
    getline(ifs, line);         //skip 1st line
    //cout << "line " << line << endl;
    while ((ifs >> t >> c >> p_val[i] >> c >> q_val[i]) && (c==',')) { i++; }


    DVD_Bundling test(2000, .03, .005, 50, 3,  0.20, 501, p_val, TSIZE); //test is instance of DVD_Bundling clas
   
    //test.print_price_array(); 
   
    test.simulate_path(TSIZE);						//run the model

    //output results in csv file
    double *revenue = test.get_revenue();				
    ofs << "t" << "," << "revenue" << endl;    
    for (t=0; t<TSIZE; t++) {
        ofs << t << "," << revenue[t] << endl;    
    }

	
    // close open files
    ifs.close();
    ofs.close();
	
    test.NPV();
	
    cout << "\nProgram ends" << endl;
}
