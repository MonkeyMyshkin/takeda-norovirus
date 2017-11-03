#include <Rcpp.h>
using namespace Rcpp;


// Helper functions for 2,3 matrix indexing
inline int index2(int i, int a, int y)
{
  return y*i + a;
}

inline void deref2(int l, int y, int & i, int & a)
{
  a = l % y;
  i = (l - a)/y;
}
// [[Rcpp::export]]
List NoroDyn(double t, NumericVector x, NumericVector params, NumericMatrix Cm, NumericVector mu, double max_age, double B, double theta)
{

  NumericVector dxdt(x.length());
  
  double alpha = params[0]; //loss of maternal antibodies
  double q = params[1]; // transmissibility
  double omega = params[2]; // seasonal amplitude
  double nu = params[3]; // scaling of asymptomatic infectiousness
  double delta = params[4]; // loss of immunity
  double epsilon = params[5]; // rate of latency loss
  double sigma = params[6]; // proportion symptomatic
  double psi = params[7]; // rate infected individuals become asymptomatic
  double gamma = params[8]; // rate asymptomatic individuals become immune

  double pi = M_PI;

  double age = 1.0/365.0; // Daily Ageing
  // (Fixed) seasonal offset for age groups
  // int ageGroupSize[7] = {7,10,8,11,13,22,10};
	
  double offset[81] = {-5.2648,-5.2648 ,-5.2648 ,-5.2648 ,-5.2648 ,-5.2648 ,-5.2648,   
					-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,
					-5.3683,-5.3683,-5.4892,-5.4892,-5.4892,-5.4892,-5.4892,-5.4892,
					-5.4892,-5.4892,-5.4201,-5.4201,-5.4201,-5.4201,-5.4201,-5.4201,
					-5.4201,-5.4201,-5.4201,-5.4201,-5.4201,-5.5237,-5.5237,-5.5237,
					-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,
					-5.5237,-5.5237,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,
					-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,
					-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,   
					-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,
					-5.6790,-5.6790};
  
  enum  state_variables {M,S,E1,E2,I,A,R,N,C};
  
  NumericVector FOI(max_age);
        
  // Calculate force of infection
  for(int a=0 ; a<max_age ; a++ )
  {
  	FOI[a] = 0.0;
    for( int a2=0; a2<max_age; a2++ )
    	{
    		double Z = (1+omega*cos(offset[a] + 2*pi*t/365.0));
        	FOI[a] += q*Z*Cm(a,a2)*(x[index2(I,a2,max_age)] + nu * x[index2(A,a2,max_age)]);
        }
  }

		dxdt[ index2(M , 0, max_age )] = B - (alpha + mu[0]) * x[index2(M,0,max_age)] - age*x[index2(M,0,max_age)];
        
        dxdt[ index2(S , 0, max_age )] = alpha*x[index2(M,0,max_age)] - (FOI[0]+ mu[0])* x[index2(S,0,max_age)] 
        				+ delta * x[index2(R,0,max_age)] - age*x[index2(S,0,max_age)];
        				
        dxdt[ index2(E1 , 0, max_age )] = FOI[0]*x[index2(S,0,max_age)] - (epsilon+mu[0]) * x[index2(E1,0,max_age)] - age * x[index2(E1,0,max_age)];
        
        dxdt[ index2(E2 , 0, max_age )] = epsilon*x[index2(E1,0,max_age)] - (epsilon+mu[0]) * x[index2(E2,0,max_age)] - age*x[index2(E2,0,max_age)];
        
        dxdt[ index2(I , 0, max_age)] = sigma*epsilon*x[index2(E2,0,max_age)] - (psi+mu[0])*x[index2(I,0,max_age)] - age*x[index2(I,0,max_age)];

        dxdt[ index2(A , 0, max_age )] = (1-sigma)*epsilon*x[index2(E2,0,max_age)] + psi*x[index2(I,0,max_age)] - (psi*gamma+mu[0]) * x[index2(A,0,max_age)]
        				 - age * x[index2(A,0,max_age)] + theta*FOI[0]*x[index2(R,0,max_age)];
        		 
        dxdt[ index2(R , 0, max_age )] = psi*gamma*x[index2(A,0,max_age)] - (delta+mu[0])*x[index2(R,0,max_age)] 
        				- age*x[index2(R,0,max_age)] - theta*FOI[0]*x[index2(R,0,max_age)];     
        
        dxdt[ index2(C , 0, max_age )] = (psi)*x[index2(I,0,max_age)];
        
        dxdt[ index2(N, 0, max_age)] = B - (age + mu[0])*(x[index2(M,0,max_age)] + x[index2(S,0,max_age)] + x[index2(E1,0,max_age)] 
        						+ x[index2(E2,0,max_age)] + x[index2(I,0,max_age)] + x[index2(A,0,max_age)] + x[index2(R,0,max_age)]);  
        
	    // Annual age groups
        for( size_t a=1 ; a<(max_age-1) ; ++a )
        {

        dxdt[ index2(M , a, max_age )] = age*(x[index2(M,a-1,max_age)] - x[index2(M,a,max_age)]) - (alpha + mu[a]) * x[index2(M,a,max_age)] ;
      
        dxdt[ index2(S , a, max_age )] = alpha*x[index2(M,a,max_age)] - (FOI[a]+ mu[a]) * x[index2(S,a,max_age)] 
        				+ delta * x[index2(R,a,max_age)] + age*(x[index2(S,a-1,max_age)] - x[index2(S,a,max_age)]);
        				
        dxdt[ index2(E1 , a, max_age )] = FOI[a]*x[index2(S,a,max_age)] - (epsilon+mu[a]) * x[index2(E1,a,max_age)] 
        				+ age * (x[index2(E1,a-1,max_age)] - x[index2(E1,a,max_age)]);
        				
        dxdt[ index2(E2 , a, max_age )] = epsilon*x[index2(E1,a,max_age)] - (epsilon+mu[a]) * x[index2(E2,a,max_age)] 
        				+ age * (x[index2(E2,a-1,max_age)] - x[index2(E2,a,max_age)]);
        				
        dxdt[ index2(I , a, max_age )] = sigma*epsilon*x[index2(E2,a,max_age)] - (psi+mu[a])*x[index2(I,a,max_age)] 
        				+ age * (x[index2(I,a-1,max_age)] - x[index2(I,a,max_age)]);
        				
        dxdt[ index2(A , a, max_age )] = (1-sigma)*epsilon*x[index2(E2,a,max_age)] + psi*x[index2(I,a,max_age)] 
        				- (psi*gamma+mu[a]) * x[index2(A,a,max_age)] + age * (x[index2(A,a-1,max_age)] - x[index2(A,a,max_age)])
        				+ theta*FOI[a]*x[index2(R,a,max_age)];

        dxdt[ index2(R , a, max_age )] = psi*gamma*x[index2(A,a,max_age)] - (delta+mu[a])*x[index2(R,a,max_age)] 
        				+ age * (x[index2(R,a-1,max_age)] - x[index2(R,a,max_age)]) - theta*FOI[a]*x[index2(R,a,max_age)];
        
        dxdt[ index2(C , a, max_age )] = (psi)*x[index2(I,a,max_age)];
           
        dxdt[ index2(N, a, max_age)] = age*(x[index2(M,a-1,max_age)] + x[index2(S,a-1,max_age)] + x[index2(E1,a-1,max_age)] + x[index2(E2,a-1,max_age)] 
        						+ x[index2(I,a-1,max_age)] + x[index2(A,a-1,max_age)] + x[index2(R,a-1,max_age)]) 
        				- (age + mu[a])*(x[index2(M,a,max_age)] + x[index2(S,a,max_age)] + x[index2(E1,a,max_age)] + x[index2(E2,a,max_age)] 
        						+ x[index2(I,a,max_age)] + x[index2(A,a,max_age)] + x[index2(R,a,max_age)]);  
       
        												
        }
        
        // Final Age Group
        dxdt[ index2(M , (max_age-1),max_age )] = age*(x[index2(M,(max_age-2),max_age)]) 
        						- (alpha +mu[max_age-1]) * x[index2(M,(max_age-1),max_age)] ;
       
        dxdt[ index2(S , (max_age-1),max_age )] = alpha*x[index2(M,(max_age-1),max_age)] 
        			- (FOI[(max_age-1)] + mu[max_age-1]) * x[index2(S,(max_age-1),max_age)] 
        				+ delta * x[index2(R,(max_age-1),max_age)] + age *(x[index2(S,(max_age-2),max_age)]);
        				
        dxdt[ index2(E1 , (max_age-1),max_age )] = FOI[(max_age-1)]*x[index2(S,(max_age-1),max_age)] 
        						- (epsilon+mu[max_age-1]) * x[index2(E1,(max_age-1),max_age)] 
        						+ age * (x[index2(E1,(max_age-2),max_age)]);
        						
        dxdt[ index2(E2 , (max_age-1),max_age )] = epsilon*x[index2(E1,(max_age-1),max_age)]
        										- (epsilon+mu[max_age-1]) * x[index2(E2,(max_age-1),max_age)] 
        										 + age * (x[index2(E2,(max_age-2),max_age)]);
        							
        dxdt[ index2(I , (max_age-1),max_age )] = sigma*epsilon*x[index2(E2,(max_age-1),max_age)] 
        						- (psi+mu[max_age-1])*x[index2(I,(max_age-1),max_age)] 
        						+ age * (x[index2(I,(max_age-2),max_age)]);
        						
        dxdt[ index2(A , (max_age-1),max_age )] = (1-sigma)*epsilon*x[index2(E2,(max_age-1),max_age)] 
        						+ psi*x[index2(I,(max_age-1),max_age)] 
        						- (psi*gamma+mu[max_age-1]) * x[index2(A,(max_age-1),max_age)] 
        						+ age * (x[index2(A,(max_age-2),max_age)]) 
        						+ theta * FOI[max_age-1] * x[index2(R,max_age-1,max_age)];
        						
        dxdt[ index2(R , (max_age-1),max_age )] = psi*gamma*x[index2(A,(max_age-1),max_age)] 
        							-(delta+mu[max_age-1])*x[index2(R,(max_age-1),max_age)] 
        							+ age * (x[index2(R,(max_age-2),max_age)]) 
        							- theta*FOI[max_age-1]*x[index2(R,max_age-1,max_age)];
        						
        dxdt[ index2(C , (max_age-1),max_age )] = (psi)*x[index2(I,(max_age-1),max_age)];	
        					   
        dxdt[ index2(N, (max_age-1),max_age)] = age*(x[index2(M,(max_age-2),max_age)] + x[index2(S,(max_age-2),max_age)] 
        								+ x[index2(E1,(max_age-2),max_age)] + x[index2(E2,(max_age-2),max_age)] 
        								+ x[index2(I,(max_age-2),max_age)] + x[index2(A,(max_age-2),max_age)] 
        								+ x[index2(R,(max_age-2),max_age)]) - (mu[max_age-1]*(x[index2(M,(max_age-1),max_age)] 
        								+ x[index2(S,(max_age-1),max_age)] + x[index2(E1,(max_age-1),max_age)] 
        								+ x[index2(E2,(max_age-1),max_age)] + x[index2(I,(max_age-1),max_age)] 
        								+ x[index2(A,(max_age-1),max_age)] + x[index2(R,(max_age-1),max_age)]));  
    
  return List::create(dxdt);
  }
  
// [[Rcpp::export]]
List NoroFP(double t, NumericVector x, NumericVector params, NumericMatrix Cm, NumericVector mu, double max_age, double B, double theta)
{

  NumericVector dxdt(x.length());
  
  double alpha = params[0]; //loss of maternal antibodies
  double q = params[1]; // transmissibility
  double omega = params[2]; // seasonal amplitude
  double nu = params[3]; // scaling of asymptomatic infectiousness
  double delta = params[4]; // loss of immunity
  double epsilon = params[5]; // rate of latency loss
  double sigma = params[6]; // proportion symptomatic
  double psi = params[7]; // rate infected individuals become asymptomatic
  double gamma = params[8]; // rate asymptomatic individuals become immune

  double pi = M_PI;

  double age = 1.0/365.0; // Daily Ageing
  // (Fixed) seasonal offset for age groups
  // int ageGroupSize[7] = {7,10,8,11,13,22,10};
	
  double offset[81] = {-5.2648,-5.2648 ,-5.2648 ,-5.2648 ,-5.2648 ,-5.2648 ,-5.2648,   
					-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,-5.3683,
					-5.3683,-5.3683,-5.4892,-5.4892,-5.4892,-5.4892,-5.4892,-5.4892,
					-5.4892,-5.4892,-5.4201,-5.4201,-5.4201,-5.4201,-5.4201,-5.4201,
					-5.4201,-5.4201,-5.4201,-5.4201,-5.4201,-5.5237,-5.5237,-5.5237,
					-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,-5.5237,
					-5.5237,-5.5237,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,
					-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,
					-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,-5.5927,   
					-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,-5.6790,
					-5.6790,-5.6790};
  
  enum  state_variables {M,S,E1,E2,I,A,R};
  
  NumericVector FOI(max_age);
        
  // Calculate force of infection
  for(int a=0 ; a<max_age ; a++ )
  {
  	FOI[a] = 0.0;
    for( int a2=0; a2<max_age; a2++ )
    	{
    		double Z = (1+omega*sin(offset[a] + 2*pi*t/365.0));
        	FOI[a] += q*Z*Cm(a,a2)*(x[index2(I,a2,max_age)] + nu * x[index2(A,a2,max_age)]);
        }
  }

		dxdt[ index2(M , 0, max_age )] = B - (alpha + mu[0]) * x[index2(M,0,max_age)] - age*x[index2(M,0,max_age)];
        
        dxdt[ index2(S , 0, max_age )] = alpha*x[index2(M,0,max_age)] - (FOI[0]+ mu[0])* x[index2(S,0,max_age)] 
        				+ delta * x[index2(R,0,max_age)] - age*x[index2(S,0,max_age)];
        				
        dxdt[ index2(E1 , 0, max_age )] = FOI[0]*x[index2(S,0,max_age)] - (epsilon+mu[0]) * x[index2(E1,0,max_age)] - age * x[index2(E1,0,max_age)];
        
        dxdt[ index2(E2 , 0, max_age )] = epsilon*x[index2(E1,0,max_age)] - (epsilon+mu[0]) * x[index2(E2,0,max_age)] - age*x[index2(E2,0,max_age)];
        
        dxdt[ index2(I , 0, max_age)] = sigma*epsilon*x[index2(E2,0,max_age)] - (psi+mu[0])*x[index2(I,0,max_age)] - age*x[index2(I,0,max_age)];

        dxdt[ index2(A , 0, max_age )] = (1-sigma)*epsilon*x[index2(E2,0,max_age)] + psi*x[index2(I,0,max_age)] - (psi*gamma+mu[0]) * x[index2(A,0,max_age)]
        				 - age * x[index2(A,0,max_age)] + theta*FOI[0]*x[index2(R,0,max_age)];
        		 
        dxdt[ index2(R , 0, max_age )] = psi*gamma*x[index2(A,0,max_age)] - (delta+mu[0])*x[index2(R,0,max_age)] 
        				- age*x[index2(R,0,max_age)] - theta*FOI[0]*x[index2(R,0,max_age)];     
        
	    // Annual age groups
        for( size_t a=1 ; a<(max_age-1) ; ++a )
        {

        dxdt[ index2(M , a, max_age )] = age*(x[index2(M,a-1,max_age)] - x[index2(M,a,max_age)]) - (alpha + mu[a]) * x[index2(M,a,max_age)] ;
      
        dxdt[ index2(S , a, max_age )] = alpha*x[index2(M,a,max_age)] - (FOI[a]+ mu[a]) * x[index2(S,a,max_age)] 
        				+ delta * x[index2(R,a,max_age)] + age*(x[index2(S,a-1,max_age)] - x[index2(S,a,max_age)]);
        				
        dxdt[ index2(E1 , a, max_age )] = FOI[a]*x[index2(S,a,max_age)] - (epsilon+mu[a]) * x[index2(E1,a,max_age)] 
        				+ age * (x[index2(E1,a-1,max_age)] - x[index2(E1,a,max_age)]);
        				
        dxdt[ index2(E2 , a, max_age )] = epsilon*x[index2(E1,a,max_age)] - (epsilon+mu[a]) * x[index2(E2,a,max_age)] 
        				+ age * (x[index2(E2,a-1,max_age)] - x[index2(E2,a,max_age)]);
        				
        dxdt[ index2(I , a, max_age )] = sigma*epsilon*x[index2(E2,a,max_age)] - (psi+mu[a])*x[index2(I,a,max_age)] 
        				+ age * (x[index2(I,a-1,max_age)] - x[index2(I,a,max_age)]);
        				
        dxdt[ index2(A , a, max_age )] = (1-sigma)*epsilon*x[index2(E2,a,max_age)] + psi*x[index2(I,a,max_age)] 
        				- (psi*gamma+mu[a]) * x[index2(A,a,max_age)] + age * (x[index2(A,a-1,max_age)] - x[index2(A,a,max_age)])
        				+ theta*FOI[a]*x[index2(R,a,max_age)];

        dxdt[ index2(R , a, max_age )] = psi*gamma*x[index2(A,a,max_age)] - (delta+mu[a])*x[index2(R,a,max_age)] 
        				+ age * (x[index2(R,a-1,max_age)] - x[index2(R,a,max_age)]) - theta*FOI[a]*x[index2(R,a,max_age)];
											
        }
        
        // Final Age Group
        dxdt[ index2(M , (max_age-1),max_age )] = age*(x[index2(M,(max_age-2),max_age)]) 
        						- (alpha +mu[max_age-1]) * x[index2(M,(max_age-1),max_age)] ;
       
        dxdt[ index2(S , (max_age-1),max_age )] = alpha*x[index2(M,(max_age-1),max_age)] 
        			- (FOI[(max_age-1)] + mu[max_age-1]) * x[index2(S,(max_age-1),max_age)] 
        				+ delta * x[index2(R,(max_age-1),max_age)] + age *(x[index2(S,(max_age-2),max_age)]);
        				
        dxdt[ index2(E1 , (max_age-1),max_age )] = FOI[(max_age-1)]*x[index2(S,(max_age-1),max_age)] 
        						- (epsilon+mu[max_age-1]) * x[index2(E1,(max_age-1),max_age)] 
        						+ age * (x[index2(E1,(max_age-2),max_age)]);
        						
        dxdt[ index2(E2 , (max_age-1),max_age )] = epsilon*x[index2(E1,(max_age-1),max_age)]
        										- (epsilon+mu[max_age-1]) * x[index2(E2,(max_age-1),max_age)] 
        										 + age * (x[index2(E2,(max_age-2),max_age)]);
        							
        dxdt[ index2(I , (max_age-1),max_age )] = sigma*epsilon*x[index2(E2,(max_age-1),max_age)] 
        						- (psi+mu[max_age-1])*x[index2(I,(max_age-1),max_age)] 
        						+ age * (x[index2(I,(max_age-2),max_age)]);
        						
        dxdt[ index2(A , (max_age-1),max_age )] = (1-sigma)*epsilon*x[index2(E2,(max_age-1),max_age)] 
        						+ psi*x[index2(I,(max_age-1),max_age)] 
        						- (psi*gamma+mu[max_age-1]) * x[index2(A,(max_age-1),max_age)] 
        						+ age * (x[index2(A,(max_age-2),max_age)]) 
        						+ theta * FOI[max_age-1] * x[index2(R,max_age-1,max_age)];
        						
        dxdt[ index2(R , (max_age-1),max_age )] = psi*gamma*x[index2(A,(max_age-1),max_age)] 
        							-(delta+mu[max_age-1])*x[index2(R,(max_age-1),max_age)] 
        							+ age * (x[index2(R,(max_age-2),max_age)]) 
        							- theta*FOI[max_age-1]*x[index2(R,max_age-1,max_age)];
    
  return List::create(dxdt);
  }
  