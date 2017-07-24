using namespace std;
#include <iostream> 
#include <fstream>

#define GIDX(i,j) (i*num_boxes+j)
#define NEW_RELEASE_WEEK 156
#define INFTY 10000
#define WEEKS_TO_BUNDLING 52
#define BUNDLING_WEEK NEW_RELEASE_WEEK + WEEKS_TO_BUNDLING 
#define TSIZE 260
class G_Array {
//private data
	private:

	public:
		double eta, phi, nu, psi, tau, step_size;
		double f;
		double *g;
		int num_boxes, prices_size;
		bool bundle;
		static const int prn_log = 1;

//public methods
	public:
		G_Array(double eta, double phi, double  nu, double psi, double tau, int num_boxes, int prices_size);
		virtual ~G_Array();
		double update_cell_purchases(int i, int j);
		double update_cell_no_purchases(int i, int j);
		double update_zero_box(int i, int j);
		virtual void initialize_array() = 0;
		double sales_revenues();
		double buy_bundle(double p_b, int i, int j);
		double buy_1(double p_t_1, int i, int j);	
		double buy_2(double p_t_2, int i, int j);
		double buy_both(double p_t_1, double p_t_2, int i, int j);
		void buy_none(int i, int j);
		
		
}; //close class definition
	
G_Array::G_Array( double eta, double phi, double  nu, double psi, double tau, int num_boxes, int prices_size): 
	eta{eta}, phi{phi}, nu{nu}, psi{psi}, tau{tau}, num_boxes{num_boxes}, prices_size{prices_size}, bundle{false} {
		step_size = psi / (num_boxes - 1);
		f = 0; 					//should be re-initialized in dervived classes
		g = new double [num_boxes * num_boxes];
		}
G_Array::~G_Array(){
	delete[] g; 
}
double inline G_Array::update_cell_purchases(int i, int j) {
	return ( (1 - phi) *(1 - nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size );
}
double inline G_Array:: update_cell_no_purchases(int i, int j) {
	return ( (1- nu) * g[GIDX(i,j)] + nu * eta * f * step_size * step_size ); 
}
double inline G_Array:: update_zero_box(int i, int j){
	return ( phi * (1 - nu) * g[GIDX(i,j)] );
}
double inline G_Array::buy_bundle(double p_b, int i, int j){
	double g_new = update_cell_purchases(i, j);
	//g[GIDX(0,0)] = update_zero_box(i, j);
	g[GIDX(0,0)] +=  phi * g[GIDX(i,j)]; 		// !!!!!!!!!!!!!!!!!!!!!!!!!!!! FIX TO INCLUDE (1-nu)
	g[GIDX(i,j)] = g_new;
	return (p_b * phi * g[GIDX(i,j)]);
}
double inline G_Array::buy_1(double p_t_1, int i, int j){
	double g_new = update_cell_purchases(i, j);
        g[GIDX(0,j)] += update_zero_box(i, j);
        g[GIDX(i,j)] = g_new;
        return (p_t_1 * phi * g[GIDX(i,j)]);
}
double inline G_Array::buy_2(double p_t_2, int i, int j){
	double g_new = update_cell_purchases(i, j);
	g[GIDX(i,0)] += update_zero_box(i, j);
	g[GIDX(i,j)] = g_new;
	return (p_t_2 * phi * g[GIDX(i,j)]);
}
double inline G_Array::buy_both(double p_t_1, double p_t_2, int i, int j){
	double g_new = update_cell_purchases(i, j);
	g[GIDX(0,0)] += update_zero_box(i, j);
	g[GIDX(i,j)] = g_new;
	return ( (p_t_1 + p_t_2) * phi * g[GIDX(i,j)] );
}
void inline G_Array:: buy_none(int i, int j){
	g[GIDX(i,j)] = update_cell_no_purchases(i, j);
	return;
}


/* -------------------------------------------------------------------------------------------------------------- */


class Bundle_Simulation {
	public:
		void simulate_path(int time);
		double compute_val_both(int i, int j);
		double compute_val_1(int i, int j);
		double compute_val_2(int i, int j);
		bool check_for_bundle(double p_b, int i, int j);
		bool bundle_check_p_1(double p_t_1, int i, int j);
		bool bundle_check_p_2(double p_t_2, int i, int j);
		bool check_1(double p_t_1, int i, int j);
		bool check_2(double p_t_2, int i, int j);
		void g_update(double p_t_1, double p_t_2,/* double  a_t_1, double a_t_2,*/ int t);
		double discount(int t);							// TBD
		double *revenue, *rev_bundle, *rev_sep;
		double omega;
		int tau;
	       	double *prices;	
		Bundle_Simulation(double omega, int tau, int prices_size);
		Bundle_Simulation(double eta, double phi, double  nu, double psi, int tau, double omega, int num_boxes, int prices_size, double* prices);
		
			

	protected:
		G_Array *g;
		double a_t_1, a_t_2;					//discount factor
		int prices_size;
};
/**/
//Bundle_Simulation::Bundle_Simulation(double omega, int tau, int prices_size): omega{omega}, tau{tau}, prices_size{prices_size} {}

/**/
Bundle_Simulation::Bundle_Simulation(double eta, double phi, double nu, double psi, int tau, double omega, int num_boxes, int prices_size, double *prices):
  	omega{omega}, tau{tau}, prices_size{prices_size}, prices{prices}, a_t_1{0}, a_t_2{0} {}
	


/**/
double Bundle_Simulation::discount (int t) {
    double a_t;
    if (t <= tau)
        a_t = 1;
    else
       a_t = omega + (1-omega) * ( (double)tau/(double)t );
    return a_t;
}


void Bundle_Simulation::simulate_path(int time){
    double p_t_1, p_t_2;                  //weekly price and dicount factor for DVD1 and DVD2
	int t;
    for (t=0; t<NEW_RELEASE_WEEK; t++) {                        //Only first DVD has been released
        p_t_1 = prices[t];
        p_t_2 = INFTY;                                          //this allows writing all function for 2 DVD case
        a_t_1 = discount(t);
        a_t_2 = 1;                                              //arbitrary placeholder value
        g_update(p_t_1, p_t_2 /*, a_t_1, a_t_2*/, t);        //update g and store the revenue for the week
    }
    for (; t < BUNDLING_WEEK; t++) {                            //new DVD released. NO BUNDLING YET
        p_t_1 = prices[t];
        //p_t_2 = INFTY;
        p_t_2 = prices[t - NEW_RELEASE_WEEK];
        a_t_1 = discount(t);
        a_t_2 = discount(t - NEW_RELEASE_WEEK);
        g_update(p_t_1, p_t_2 /*, a_t_1, a_t_2*/, t);
    }
    for (; t < time; t++){                                      //Bundle offered. Continue until end of tracking period
        g->bundle = true;          //bundle se calculation will now be activated
        p_t_1 = prices[t];
        //p_t_2 = INFTY;
        p_t_2 = prices[t - NEW_RELEASE_WEEK];
        a_t_1 = discount(t);
        a_t_2 = discount(t - NEW_RELEASE_WEEK);
        g_update(p_t_1, p_t_2 /*, a_t_1, a_t_2*/, t);
    }
    return ;
}

double inline Bundle_Simulation::compute_val_both(int i, int j){
	return (a_t_1*((double)i + 0.5) + a_t_2 *((double) j + 0.5))* g->step_size;
}
double inline Bundle_Simulation::compute_val_1(int i, int j){
	return a_t_1 *((double) i + 0.5) * g->step_size;
}
double inline Bundle_Simulation::compute_val_2(int i, int j){
	return a_t_2 *((double) j + 0.5) * g->step_size;
}
bool inline Bundle_Simulation::check_for_bundle(double p_b, int i, int j) {
	double val_both = compute_val_both(i, j);
	//return less_equal<double>(p_b, val_both);						// TBD
	return (p_b <= val_both ? 1 : 0); 
}
bool inline Bundle_Simulation::bundle_check_p_1(double p_t_1, int i, int j){
	double val_1 = compute_val_1(i, j);
	//return less(val_1, (p_t_1 - 5) );
	return (val_1 < (p_t_1 - 5) ? 1 : 0);
}
bool inline Bundle_Simulation::bundle_check_p_2(double p_t_2, int i, int j){
	double val_2 = compute_val_2(i, j);
	//return less(val_2, (p_t_2 - 5) );
	return (val_2 < (p_t_2 - 5) ? 1:0 );
}
bool inline Bundle_Simulation::check_1(double p_t_1, int i, int j){
	double val_1 = compute_val_1(i, j);
	//return less_equal (val_1, p_t_1);
	return (p_t_1 <= val_1 ? 1:0);
}	
bool inline Bundle_Simulation::check_2(double p_t_2, int i, int j){
	double val_2 = compute_val_2(i, j);
	//return less_equal(val_2, p_t_2);
	return (p_t_2 <=  val_2 ? 1:0);
}

void Bundle_Simulation::g_update(double p_t_1, double p_t_2 /*, double  a_t_1, double a_t_2*/, int t){
	
	double temp = 0;
	double rev =0;
	double p_b = p_t_1 + p_t_2 - 5;
	for (int i = 0; i < g->num_boxes; i++){
		for (int j = 0; j < g->num_boxes; j++) {
			if ( check_for_bundle(p_b, i, j) && ( g->bundle == true ) ) {
				if ( bundle_check_p_1(p_t_1, i, j) ) {
					temp = g->buy_2(p_t_2, i, j);
					rev += temp;
				}
				else if ( bundle_check_p_2(p_t_2, i, j) ) {
					temp = g->buy_1(p_t_1, i, j);
					rev += temp;
				}	
				else {
					rev += g->buy_bundle(p_b, i, j);
				}
			}
			else if ( check_1(p_t_1, i, j) ) {
				if ( check_2(p_t_2, i, j) ) {
					temp = g->buy_both(p_t_1, p_t_2, i, j);
					rev += temp;
				}
				else { 
					temp = g->buy_1(p_t_1, i, j);
					rev += temp;

				}
			}
			else if ( check_2(p_t_2, i, j) ) {
				temp = g->buy_2(p_t_2, i, j);
				rev += temp;
			}
			else{ 
				g->buy_none(i, j);
			}
		} //close j loop
	} //close i loop
	revenue[t] = rev;
}
	

class Independent:public G_Array {
	public:
		Independent(double eta, double phi, double  nu, double psi, double tau, int num_boxes, int prices_size);
	//	~Independent(){;};
		void initialize_array();
		;
		
		

};

Independent::Independent(double eta, double phi, double  nu, double psi, double tau, int num_boxes, int prices_size):
	 G_Array(eta, phi, nu, psi, tau, num_boxes, prices_size){
	f = (double)1/(psi*psi);
	initialize_array();
	
}
//Independent::~Independent():~G_Array(){;}
void Independent::initialize_array(){  

int i, j;
for (i=0; i < num_boxes; i++) { 
	for (j=0; j < num_boxes; j++) {
                 if( (i == num_boxes -1) || (j == num_boxes -1) ) {
                         g[GIDX(i,j)] = 0;
			//cout<< "a" << endl;
                        }
                 else{
                         g[GIDX(i,j)] = eta * f * step_size*step_size;
                        }                       // number of people in each box of the approxmation for the g(v)
        }
    }


}


/**/
class Independent_Simulation: public Bundle_Simulation {
	public:
		Independent_Simulation(double eta, double phi, double nu, double psi, int tau, double omega, int num_boxes, int prices_size, double *prices);
		~Independent_Simulation(){}
		void NPV();
};

/**/
Independent_Simulation::Independent_Simulation(double eta, double phi, double nu, double psi, int tau, double omega, int num_boxes, int prices_size, double *prices):
	Bundle_Simulation(eta, phi, nu, psi, tau, omega, num_boxes, prices_size, prices) {
		revenue = new double[prices_size];
		//rev_bundle = new double[prices_size];
		//rev_sep = new double[prices_size];
	        g = new Independent(eta, phi, nu, psi, tau,  num_boxes, prices_size);
}
	
/**/
void Independent_Simulation::NPV(){
	double total = 0;
	for (int i = 0; i<TSIZE; i++){
		total += revenue[i];
		cout<<revenue[i] <<endl;
	}
	cout<<total<<endl;
	return;
}

/**/
int main() {
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

    ofs.open("weekly_revenues.csv");
    if (! ofs.is_open() ) {
        cout << "Error. Cannot open output file" << endl;
        return 4;
    }

    // read input data
    getline(ifs, line);         //skip 1st line
    //cout << "line " << line << endl;
    while ((ifs >> t >> c >> p_val[i] >> c >> q_val[i]) && (c==',')) { i++; }



	
	/* Independent * test;
	test = new Independent (2000, .03, .005, 50, 3, 10, 260);
	test->initialize_array(); 
	delete test; */

	Independent_Simulation *test2 = new Independent_Simulation(2000, .03, .005, 50, 3,  0.20, 501, TSIZE, p_val);
	test2->simulate_path(260);
	test2->NPV();
	delete test2;
	
	ifs.close();
	ofs.close();


	return 0;
}
 
		
		

		
		


	
			
			
	
