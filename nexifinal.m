
data = readtable('C:\Users\USER\Downloads\NEXI.MI.csv');

dates = datetime(data.Date, 'InputFormat', 'yyyy-MM-dd');
prices = data.AdjClose;

figure;
plot(dates, prices);
title('NEXI.MI Historical Prices');
xlabel('Date');
ylabel('Adjusted Close Price (€)');
grid on;

returns = diff(log(prices));  

volatility = std(returns) * sqrt(252);
fprintf('Annualized Volatility: %.2f%%\n', volatility * 100);

n_simulations = 50000;
T = 1;  
n_steps = 252;  
dt = T/n_steps;
S0 = prices(end);  

simulated_prices = zeros(n_simulations, n_steps);
for i = 1:n_simulations
    simulated_prices(i, 1) = S0;
    for j = 2:n_steps
        dW = sqrt(dt) * randn;
        simulated_prices(i, j) = simulated_prices(i, j-1) * exp((0.04 - 0.5 * volatility^2) * dt + volatility * dW);
    end
end


figure;
plot(1:n_steps, simulated_prices(1:10, :));
title('Simulated Price Paths for NEXI.MI');
xlabel('Time Steps');
ylabel('Price (€)');
grid on;


premium = 0.765;  
contracts = 3;    
lot_size = 100;   
K = 6.00;         % Strike Price
T = (datenum('15-Mar-2025') - now) / 365;  % Time to maturity in years
sigma = volatility;  
S0 = prices(end);  
risk_free_rates = linspace(0.03, 0.039, 5);  % Risk-free rate scenarios from 3% to 3.9%


total_cost = premium * contracts * lot_size;


option_prices = zeros(1, length(risk_free_rates));


for j = 1:length(risk_free_rates)
    r = risk_free_rates(j);
    
    d1 = (log(S0/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T));
    d2 = d1 - sigma*sqrt(T);
    
    % Calculate the call option price using Black-Scholes formula
    call_price = S0*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
    option_prices(j) = call_price;
end


disp('Option Prices under Different Risk-Free Rates:');
disp(table(risk_free_rates', option_prices', 'VariableNames', {'Risk_Free_Rate', 'Option_Price'}));


max_profit = inf;  
max_loss = total_cost;  
break_even_price = K + premium;  


fprintf('Total Cost: €%.2f\n', total_cost);
fprintf('Maximum Loss: €%.2f\n', max_loss);
fprintf('Break-even Stock Price at Expiration: €%.2f\n', break_even_price);




stock_prices = linspace(K * 0.5, K * 2, 100);


PnL_euros = max(stock_prices - K, 0) * contracts * lot_size - total_cost;


PnL_percentage = (PnL_euros / total_cost) * 100;


figure;
plot(stock_prices, PnL_percentage, 'b-', 'LineWidth', 2);
hold on;
yline(0, 'k--');  
xline(break_even_price, 'r--', 'Break-even');
title('PnL in Percentage for 6.00 Call Option (Maturity: March 2025)');
xlabel('Stock Price at Expiration (€)');
ylabel('PnL (%)');
grid on;

% Plot PnL in Euros
figure;
plot(stock_prices, PnL_euros, 'g-', 'LineWidth', 2);
hold on;
yline(0, 'k--');  
xline(break_even_price, 'r--', 'Break-even');
title('PnL in Euros for 6.00 Call Option (Maturity: March 2025)');
xlabel('Stock Price at Expiration (€)');
ylabel('PnL (€)');
grid on;

fprintf('Total Cost: €%.2f\n', total_cost);
fprintf('Maximum Loss: €%.2f\n', max_loss);
fprintf('Break-even Stock Price at Expiration: €%.2f\n', break_even_price);
