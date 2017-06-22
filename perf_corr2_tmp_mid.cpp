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
	double npv;  //These vars are used to compute seperate NPVs
	double npv_sep;
	double npv_bundle;

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
    int bundling_week;  // /week at which bundling begins

    double step_size;        
    double f;            // uniform density 
    double *g;           // The number of people with a valuations in intervals corresponding to coarsness of mesh
                        // it's a variable size array
                            
    double* revenue;	// pointer to array storing weekly revenue
    double* rev_bundle; 
    double* rev_sep;
    int revenue_size;	//number of weeks in revenue array
    
    bool bundle;	// true if bundling, false if no bundling
    static const int prn_log = 0;

//public methods
public:
    DVD_Bundling(double _eta, double _phi, double _nu, double _psi, int _tau, double _omega, 
                 int _num_boxes, double *_prices, int _prices_size, int _bundling_week);
    ~DVD_Bundling();
    void simulate_path(int time);           //The main to generate model results
    double* get_revenue(){return revenue;}   // allows use of revenue array outside class
    double* get_rev_bundle(){return rev_bundle;}
    double* get_rev_sep(){return rev_sep;}
    void  comp_NPV();	 		    //Compute Net Present Value at t = 0	
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
    void g_update(double p_t_1, double p_t_2, double  a_t_1, double a_t_2, int t);	//This method computes the new g(v) and then computes weekly sales
    double discount(int t);                         				//This function computes the discount rate for each time value
};

/* constructor */
DVD_Bundling::DVD_Bundling(double _eta, double _phi, double _nu, double _psi, int _tau, double _omega, 
                           int _num_boxes, double *_prices, int _prices_size, int _bundling_week)
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
    bundle = false;
    step_size = psi/(num_boxes - 1);
    f = (double)1/(psi);      //uniform density
	
    g = new double[num_boxes * num_boxes]; 

    for (i=0; i < num_boxes; i++) {
        for (j=0; j < num_boxes; j++) {
			if( (i == num_boxes -1) || (j == num_boxes -1) ) {
				g[GIDX(i,j)] = 0;
			}
			else if( (i ==j) ){	
				g[GIDX(i,j)] = eta * f * step_size;
			}
			else{ 
				g[GIDX(i,j)] = 0;
			}			// number of people in each box of the approxmation for the g(v)
        }
    } 
   //[GIDX(num_boxes-1,num_boxes-1)] = 0;                      //right bottom element
    revenue = new double[_prices_size]; 				//allocate memory for revenue array and assign to revenue pointer
    rev_bundle = new double[_prices_size];
    rev_sep = new double[_prices_size];
    revenue_size = _prices_size; 				//size of revenue array
}

/* desctructor */
DVD_Bundling::~DVD_Bundling()
{
    if (prn_log)						//used for debugging
        cout << "DVD_Bundling::~DVD_Bundling()" << endl;

    delete[] g;							//free memory allocated for g array
    delete[] revenue;						//free memory allocated for revenue array
	delete[] rev_bundle;
	delete[] rev_sep;

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
    	g_update(p_t_1, p_t_2, a_t_1, a_t_2, t); 	//update g and store the revenue for the week
    }
    for (; t < bundling_week; t++) {				//new DVD released. NO BUNDLING YET
		p_t_1 = prices[t];
		//p_t_2 = INFTY;
		p_t_2 = prices[t - NEW_RELEASE_WEEK];	
 		a_t_1 = discount(t);
		a_t_2 = discount(t - NEW_RELEASE_WEEK);
		g_update(p_t_1, p_t_2, a_t_1, a_t_2, t);
    }
    for (; t < time; t++){					//Bundle offered. Continue until end of tracking period
		bundle = true;		//bundle se calculation will now be activated
		p_t_1 = prices[t];
		//p_t_2 = INFTY; 
    	p_t_2 = prices[t - NEW_RELEASE_WEEK];
		a_t_1 = discount(t);
		a_t_2 = discount(t - NEW_RELEASE_WEEK);
		g_update(p_t_1, p_t_2, a_t_1, a_t_2, t); 
    }
    return ;
}

/* */
double DVD_Bundling::buy_bundle(double p_b, int i, int j){
	double g_new;
	if ( (i==j) ){
		g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size;
	}
	else {
		g_new = (1-phi) * (1-nu)* g[GIDX(i,j)];
	}
	g[GIDX(0,0)] += phi * g[GIDX(i,j)];							//Already bought both
	g[GIDX(i,j)] = g_new;
	return (p_b * phi * g[GIDX(i,j)]);
}

/* */
double DVD_Bundling::buy_1 (double p_t_1, int i,int j){
	double g_new;
	if ( (i==j) ){
		g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size;
	}
	else{
		g_new = (1 - phi) * (1-nu) * g[GIDX(i,j)]; 
	}
	g[GIDX(0,j)] += phi * (1 - nu) * g[GIDX(i,j)];							//May buy DVD2  in future
	g[GIDX(i,j)] = g_new;
	return (p_t_1 * phi * g[GIDX(i,j)]);
}

double DVD_Bundling::buy_2 (double p_t_2,int i,int j){
	double g_new;
	if ( (i==j) ){
		g_new = (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size;
	}
	else{
		g_new = (1-phi) * (1-nu) * g[GIDX(i,j)];
	}
	g[GIDX(i,0)] += phi * (1 - nu) * g[GIDX(i,j)];							//May buy DVD1 in future
	g[GIDX(i,j)] = g_new;
	return  (p_t_2 * phi * g[GIDX(i,j)]);
}

double DVD_Bundling::buy_both(double p_t_1, double p_t_2, int i, int j){
	double g_new;
	if ( (i==j) ){
		g_new = (1- phi) * (1 - nu) * g[GIDX(i,j)] + nu * eta * f *step_size;
	}
	else{
		g_new = (1-phi) * (1-nu) * g[GIDX(i,j)];
	}
	g[GIDX(0,0)] += phi * (1 - nu) * g[GIDX(i,j)]; 
	g[GIDX(i,j)] = g_new;
	return ( (p_t_1 + p_t_2) * phi * g[GIDX(i,j)] );
}
void DVD_Bundling::buy_none (int i, int j){
	double g_new;
	if ( (i==j) ){ 
		g_new =(1- nu) * g[GIDX(i,j)] + nu * eta * f* step_size;
	}
	else {
		g_new = (1-nu) * g[GIDX(i,j)];
	}
	g[GIDX(i,j)] = g_new;
	return;
} 
			
/* */
void  DVD_Bundling::g_update(double p_t_1, double p_t_2, double  a_t_1, double a_t_2, int t)
{
    int i;
    int j;
	double temp = 0;
	double rev_bundling = 0;
	double rev_seperates = 0;
    double rev = 0;
    //double p_b = INFTY;
    double p_b = p_t_1 + p_t_2 - 5;										//price of bundle
    for (i = 0; i < num_boxes; i++){
	    for (j = 0; j < num_boxes; j++){
		    if ( (p_b <= (a_t_1*((double) i + 0.5) + a_t_2 *((double) j + 0.5))* step_size) && ( bundle == true) ){			//Buy bundle ?
			
				if ( (a_t_1 * ((double) i + 0.5) * step_size) < (p_t_1 - 5) ){
 					temp = buy_2(p_t_2, i, j);
					rev += temp;
					rev_seperates += temp;
				}
				else if ( (a_t_2 *((double) j + 0.5) * step_size) < (p_t_2 -5) ){
					temp = buy_1(p_t_1, i, j);
					rev += temp;
					rev_seperates += temp;
				}
				else {
					temp = buy_bundle(p_b, i, j);
					rev += temp;
					rev_bundling += temp;
				}
			}
		
		else if (p_t_1 <= a_t_1 *((double) i  + 0.5)* step_size){							//Buy only DVD1 ?
			if (p_t_2 <= a_t_2 *((double) j + 0.5) * step_size){							//Buy only DVD1 ?
				temp = buy_both(p_t_1, p_t_2, i, j);
				rev += temp;
				rev_seperates += temp;
			}
			else{
				temp = buy_1(p_t_1, i, j);
				rev += temp;
				rev_seperates += temp;
			}


		}
		else if(p_t_2 <= a_t_2  * ((double) j + 0.5) * step_size) {
			temp = buy_2(p_t_2, i, j);
			rev += temp;									//But only DVD2 ?
			rev_seperates += temp; 
		}
		

		else {												//Buy nothing.
			buy_none(i, j);
			
		} 	

	}			//close j loop
  }			//close i loop
  
	revenue[t] = rev;
	rev_bundle[t] = rev_bundling;
	rev_sep[t] = rev_seperates;	
	return; 
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
void DVD_Bundling::comp_NPV()
{
    int t; 
    double r = .02; //discount rate
    npv = 0; //Net Present Value
	npv_sep = 0;
	npv_bundle =0;
    for (t = 0; t<TSIZE; t++){
	    npv += revenue[t] * (1/pow((1+r),(double) t));
    	npv_sep += rev_sep[t] * (1/pow((1+r),(double) t));
		npv_bundle += rev_bundle[t] * (1/pow((1+r),(double) t));
	}
    //cout <<"\nNPV: " << npv << endl;
    return;
}


/* */
//#define TSIZE 260   //table size - num lines in the input file
int main()
{
    double q_val[TSIZE];
    double p_val[TSIZE];

	double npv;
	double npv_sep;
	double npv_bundle;

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

    ofs.open("perf_corr_week_iteration.csv");
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
        pTest = new DVD_Bundling(2000, .03, .005, 50, 3,  0.20, 1501, p_val, TSIZE, t);
        pTest->simulate_path(TSIZE);        //run the model
        pTest->comp_NPV();
		npv = pTest->npv;
		npv_sep = pTest->npv_sep;
		npv_bundle = pTest->npv_bundle;
        cout << "t=" << t-NEW_RELEASE_WEEK << " NPV: " << npv << endl;    
        ofs << t-NEW_RELEASE_WEEK << "," << npv << "," << npv_sep << "," << npv_bundle <<  endl;
        delete pTest; 
    }

/*
    //output results in csv file
    double *revenue = pTest->get_revenue();				
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
