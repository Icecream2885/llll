$title ADMM-based solution procedure for Urban Layout Problem
set i building /101,102,201,202/
hb(i) home building /101,102/
ob(i) office building /201,202/
;
set k location /301,302,401,402/
hl(k) /301,302/
ol(k) /401,402/
;
alias(i,j);
alias(k,l);

parameter c(k,l) travel cost between location k and l /
301.401  2,
401.301  1000,
301.402  3,
402.301  1000,
302.401  1,
401.302  1000,
302.402  4,
402.302  1000,
301.302  1000,
302.301  1000,
401.402  1000,
402.401  1000
/;

parameter CAP_IN(j)  capacity of office building /
201  7,
202  5
/;
parameter CAP_OUT(i)  capacity of home building /
101  4,
102  6
/;

set n iteration index /1*10/;

parameter Upper_bound;
Upper_bound=10000000;

parameter Lower_bound;
Lower_bound=-10000000;

parameter subgradient_k(k);
subgradient_k(k)=0;
*parameter stepsize;
*stepsize=0.5;
parameter subgradient_i(i);
subgradient_i(i)=0;
parameter subgradient_ho;
subgradient_ho=0;

parameter lamda_k(k);
*LR multiplier;
lamda_k(k)=200;
parameter lamda_i(i);
lamda_i(i)=200;
parameter lamda_ho;
lamda_ho=200;

parameter rho;
rho=2000;

parameter gap;
gap=0;

parameter value(n) /
1  1
/;

binary variables
x(i,k) whether or not building i is assigned to location k in the primal model
x_lb(i,k)   lower bound model
*x_ub(i,k)
x_location(i,k)   ADMM model of sub location problem
;

positive variable
f(i,j)  traffic flow between building i and building j in the primal model
f_lb(i,j)   lower bound model
*f_ub(i,j)
f_flow(i,j)  ADMM model of sub flow problem
;

variables
z total travel cost
z_lb lower bound
z_up upper bound
z_location  get cost when traffic flow is fixed in the admm model of sub location problem
z_flow      get cost when location desicion is fixed in the admm model of sub flow problem
;

equations
cost   define objective function of the primal problem
assign1(k)   each location k only can be placed with one building
assign2(i)   each building i only can be assigned to one location
total_travel_demand    total travel demand bewteen home building and office building is a fixed value
capacity_in(j)     building's attraction capacity constraints in the primal probelm
capacity_out(i)    building's generating capacity constraints in the primal problem
LB   define objective function of lower bound of the problem
Sub_location  define objective function of admm model of sub location problem
Sub_flow    define objective function of admm model of sub flow problem
capacity_sub_in(j)   building's attraction capacity constraints in the admm model of sub flow probelm
capacity_sub_out(i)  building's generating capacity constraints in the admm model of sub flow probelm
;

cost..  z=e=sum((i,j,k,l),f(i,j)*c(k,l)*x(i,k)*x(j,l));
assign1(k)..  sum(i,x(i,k))=e=1;
assign2(i)..  sum(k,x(i,k))=e=1;
total_travel_demand..  sum((i,j)$(hb(i) and ob(j)),f(i,j))=e=10;
capacity_in(j)$(ob(j))..  sum(i$(hb(i)),f(i,j))=l=CAP_IN(j);
capacity_out(i)$(hb(i)).. sum(j$(ob(j)),f(i,j))=l=CAP_OUT(i);
LB..   z_lb=e=sum((i,j,k,l),f_lb(i,j)*c(k,l)*x_lb(i,k)*x_lb(j,l))+sum(k$(hl(k)),lamda_k(k)*(sum(i,x_lb(i,k))-1))+sum(i,lamda_i(i)*(sum(k,x_lb(i,k))-1))+lamda_ho*(sum((i,j)$(hb(i) and ob(j)),f_lb(i,j))-1);
Sub_location..   z_location=e=sum((i,j,k,l),f_flow.l(i,j)*c(k,l)*x_location(i,k)*x_location(j,l))+sum(k,lamda_k(k)*(sum(i,x_location(i,k))-1))+sum(i,lamda_i(i)*(sum(k,x_location(i,k))-1))+rho*0.5*(sum(k,abs(sum(i,x_location(i,k))-1)**2)+sum(i,abs(sum(k,x_location(i,k))-1)**2));
Sub_flow..   z_flow=e=sum((i,j,k,l),f_flow(i,j)*c(k,l)*x_location.l(i,k)*x_location.l(j,l))+lamda_ho*(sum((i,j)$(hb(i) and ob(j)),f_flow(i,j))-10)+rho*0.5*abs(sum((i,j)$(hb(i) and ob(j)),f_flow(i,j))-10)**2;
capacity_sub_in(j)$(ob(j))..    sum(i$(hb(i)),f_flow(i,j))=l=CAP_IN(j);
capacity_sub_out(i)$(hb(i))..   sum(j$(ob(j)),f_flow(i,j))=l=CAP_OUT(i);

*standard urban layout model
model urban_layout  /cost,assign1,assign2,total_travel_demand,capacity_in,capacity_out/;
solve urban_layout using minlp minimizing z;
display x.l,f.l,z.l;

*ADMM model of urban layout
*sub_location_problem
model sub_location_model  /Sub_location,capacity_sub_in,capacity_sub_out/;
*display  x_location(i,k);

*sub_flow_problem
model sub_flow_model  /Sub_flow,capacity_sub_in,capacity_sub_out/;
*display  f_flow(i,j);

*lower bound model
model lower_bound_urban_layout /LB,capacity_in,capacity_out/;
*display x_lb(i,k),f_lb(i,j);


loop(n,
     if(value(n)=1,
*original traffic flow f_flow
        f_flow.l('101','201')=2;
        f_flow.l('101','202')=2;
        f_flow.l('102','201')=3;
        f_flow.l('102','202')=3;
        solve sub_location_model using minlp minimizing z_location;
        solve sub_flow_model using dnlp minimizing z_flow;
        subgradient_k(k)= max(0,(sum(i,x_location.l(i,k))-1));
        subgradient_i(i)= max(0,(sum(k,x_location.l(i,k))-1));
        subgradient_ho = max(0,(sum((i,j)$(hb(i) and ob(j)),f_flow.l(i,j))-10));
        lamda_k(k)= max(0,lamda_k(k)+rho*subgradient_k(k));
        lamda_i(i)= max(0,lamda_i(i)+rho*subgradient_i(i));
        lamda_ho = max(0,lamda_ho+rho*subgradient_ho);
        z_up.l = sum((i,j,k,l),f_flow.l(i,j)*c(k,l)*x_location.l(i,k)*x_location.l(j,l));
        Upper_bound = min(Upper_bound,z_up.l);
        solve lower_bound_urban_layout using minlp minimizing z_lb;
        Lower_bound = max(Lower_bound,z_lb.l);
        gap=(Upper_bound-Lower_bound)/Upper_bound;
        display Lower_bound;
        display Upper_bound;
        display gap;
     else
        solve sub_location_model using minlp minimizing z_location;
        solve sub_flow_model using dnlp minimizing z_flow;
        subgradient_k(k)= max(0,(sum(i,x_location.l(i,k))-1));
        subgradient_i(i)= max(0,(sum(k,x_location.l(i,k))-1));
        subgradient_ho = max(0,(sum((i,j)$(hb(i) and ob(j)),f_flow.l(i,j))-10));
        lamda_k(k)= max(0,lamda_k(k)+rho*subgradient_k(k));
        lamda_i(i)= max(0,lamda_i(i)+rho*subgradient_i(i));
        lamda_ho = max(0,lamda_ho+rho*subgradient_ho);
        z_up.l = sum((i,j,k,l),f_flow.l(i,j)*c(k,l)*x_location.l(i,k)*x_location.l(j,l));
        Upper_bound = min(Upper_bound,z_up.l);
        solve lower_bound_urban_layout using minlp minimizing z_lb;
        Lower_bound = max(Lower_bound,z_lb.l);
        gap=(Upper_bound-Lower_bound)/Upper_bound;
        display lamda_k;
        display lamda_i;
        display lamda_ho;
        display Lower_bound;
        display Upper_bound;
        display gap;
        display x_lb.l;
        display f_lb.l;
        display x_location.l;
        display f_flow.l;
        display z_up.l;
     );
);


File Urban_Layout_ADMM/Urban_Layout_ADMM.dat/;

put Urban_Layout_ADMM;

loop((i,k),put @1, i.tl, @10, k.tl, @30, x_location.l(i,k)/);

loop((i,j),put @1, i.tl, @10, j.tl, @30, f_flow.l(i,j)/);








