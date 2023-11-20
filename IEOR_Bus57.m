clear
clc


[Ybus, Yf, Yt] = makeYbus(case57); 
Y=full(imag(Ybus));
L=-Y;
for i=1:57
    L(i,i)=0;
    L(i,i)=-sum(L(i,:)); 
end

%% Q1


rng('default')
i = 1; j = 1;
n = 57; norm = 10000;
demand_norm = 26869/norm;

generator_index =  [1 2 3 6 8 9 12];
random_num = rand(50,1);
total = sum(random_num );
PD = zeros(n,1);
PD(setdiff(1:n, generator_index)) = (random_num/total)* demand_norm;

indexer =[1 2 3 4 5 6 7]; 
indexer_PG  =zeros(n,1);
indexer_PG(generator_index) =1;

cvx_begin quiet
    variables PG(n) theta(n)
    minimize  indexer * PG(generator_index)
    subject to

        PG .* (1 - indexer_PG) == 0;
        PG - PD == L * theta;
       
        for  i = 1:n
            for j = 1:n
                if sign(L) ~= 0
                    abs(theta(i) - theta(j)) <= pi/10;
                end
            end
        end
        abs(theta(4) -  theta(6)) <= 0;
        abs(theta(8) -  theta(9)) <= 0;
        PG >= 0;
 cvx_end
Q1_thetas = theta;
Q1_costs = indexer * PG(generator_index) * norm
Q1_PG = PG * norm



%% Q2

PG_all = zeros(n,24);
cost_all = zeros(1,24);
lambda_all = [];
demand_list = ([20477 19854 19223 18902 18973 19687 21188 22541 22070 20910 ...
    19115 17753 17151 17166 17604 18392 19667 21663 23959 26034 26869 26126 24958 23422]);
demand_norm_list = demand_list/norm;

random_num_list = rand(50,24);
total_list = sum(random_num_list );
a = zeros(1,24);
PD_50= (random_num_list/total_list)* demand_norm_list;
PD_list = [a;a;a;PD_50(1:2,:);a;PD_50(3,:);a;a;PD_50(4:5,:);a;PD_50(6:50,:)];


for hour = 1:24

    cvx_begin quiet
        variables PG_2(n) theta_2(n)
        dual variables lambda;
        minimize  indexer * PG_2(generator_index)+0.3*(PG_2(generator_index)' *PG_2(generator_index))
        subject to
    
            PG_2 .* (1 - indexer_PG) == 0;
            lambda: PG_2 - PD_list(:,hour) == L * theta_2;
           
            for  i = 1:n
                for j = 1:n
                    if sign(L) ~= 0
                        abs(theta_2(i) - theta_2(j)) <= pi/10;
                    end
                end
            end
            abs(theta_2(4) -  theta_2(6)) <= 0;
            abs(theta_2(8) -  theta_2(9)) <= 0;
            PG_2 >= 0;
     cvx_end

     PG_all(:,hour) = PG_2;
     cost_all(hour) = indexer * PG_2(generator_index)+0.3*(PG_2(generator_index)' *PG_2(generator_index));
     
     lambda_all =[lambda_all lambda];

end
Q2_PG = PG_2 * norm;
Q2_cost_all = cost_all;
Q2_lambda_all = lambda_all;
disp('The electricity prices at all buses at each of these times:')
lambda_all'

%% Q2 -- selected bus 5,23,36 and 48 and plot the prices over time
selected_buses = [5,23,36,48];
price_seleted_buses = Q2_lambda_all(selected_buses,:);
x =  1:24;
for i = 1:4   
    yi = price_seleted_buses(i,:);
    plot(x,yi)    
    hold on    
end

xlabel('Hour');
ylabel('Price');
legend('Bus 5','Bus 23','Bus 36','Bus 48')
title('Prices for selected Bus')
hold off
figure()

disp('There are two price peaks and two price troughs in each 24-hour cycle.') 
disp('The first price peak occurs roughly at 8am, and the second price peak occurs roughly at 9pm.')

disp('The first price trough occurs roughly between 4 and 5am, and the second price trough occurs roughly between 1 and 2pm.') 
disp('Each price peak is sharp, but the price troughs can be flat and last for approximately an hour. ')
disp('The average rate of price change from the afternoon to midnight is greater than that from midnight to noon. ')
disp('The increase from trough hours to peak hours is almost symmetric to the decrease from peak hours to trough hours. ')
disp('Finally, the difference in price between the second trough and the peak is more than twice that of the first trough.')

%% Q2 -- max and min price 
[max_price, max_index] = max(Q2_lambda_all(:));
[max_row,max_col] = ind2sub(size(Q2_lambda_all),max_index);
max_price_bus = max_row;

[min_price, min_index] = min(Q2_lambda_all(:));
[min_row,min_col] = ind2sub(size(Q2_lambda_all),min_index);
min_price_bus = min_row;

min_p = sprintf('The buses %d with the lowest prices %d',min_row  ,round(min_price,3));
max_p = sprintf('The buses %d with the highest prices %d',max_row  ,round(max_price,3) );
disp(min_p);
disp(max_p);


mean_prices = mean(Q2_lambda_all, 2);
min_index = find(mean_prices == min(mean_prices)); % Bus 4 -> it is connected to the edge with capacity 0
max_index = find(mean_prices == max(mean_prices)); % Bus 9 -> it is connected to the edge with capacity 0

hold on
yyaxis left
plot(demand_norm_list)
yyaxis right
plot(Q2_lambda_all(min_index, :))
grid
legend("Normalized demand", "Price at min cost bus") 
hold off
figure()
hold on
yyaxis left
plot(demand_norm_list)
yyaxis right
plot(Q2_lambda_all(max_index, :))
legend("Normalized demand", "Price at max cost bus") 
hold off


disp('The min and max price of the buses connected to an edge with capacity zero. ')
disp(['Prices of both the min and max cost bus seems to vary directly ' ...
    'proportionately with the current demand.'])
disp('Interesting to note that the price of the min cost bus is always negative.')

%%  Q3
P = zeros(n,n);

x_ij = zeros(n,n);
for i = 1:n
    for j = 1:n
        if L(i,j) ~= 0 
            x_ij(i,j) = -1/L(i,j);
        end
    end
end 

x_ij(x_ij == diag(x_ij)) == 0;

for row = 1:n
    for col = 1:n
        if L(row,col) ~= 0 
            P(row,col) = (theta(row) - theta(col)) ./ (x_ij(row,col));
        end
    end
end

bus_err_all = zeros(n,3);
edge_err_all = zeros(n,3*n);

for z = 1:3

    PG_hat = PG;
    rand_idx = randi(n,1,z);
    PG_hat(rand_idx) = PG_hat(rand_idx) + randi([101 10000],1,z)';
    injections_idx = rand_idx;

    cvx_begin quiet
        variables e(n) theta_3(n) edge_err(n,n)
        minimize (sum(abs(e)) + sum(sum(abs(edge_err))))
        subject to
    
            PG_hat + e  == sum(P,2);  % P_hat + e  = P
                  
            for  i = 1:n
                for j = 1:n
                    if sign(L) ~= 0
                        P(i,j) + edge_err(i,j) == abs(theta_3(i) -theta_3(j)) ./ x_ij;
                    end
                end
            end
           
     cvx_end
     
     round_error = round(e,0);

     disp('The attact nodes');
     rand_idx'
     disp('All buses error');
     round_error

     bus_err_all(:,z) = e;
     edge_err_all(:,z*(1:57)) =edge_err;

     if any(edge_err ~= 0)
         disp('non-zero edge error detected');
         [row_att, col_att] = find(edge_err);
     else 
         disp('All edge errors are zero');
     end
        

end











