% B�NYAM�N BERAT GEZER  210102002061 



% 1 soru a ��kk�

n = 0:50;  % ayr�k zaman aral��� 
x_n = exp(-(0.05 + 0.3j) * n);  % ayr�k zaman i�aretinin  denklemi


% a ��kk� i�in grafikler

figure;  

subplot(2,2,1);
stem(n, real(x_n));  % Ger�ek k�s�m i�in 
title('Ger�ek K�s�m');
xlabel('n');
ylabel('Re{x[n]}');

subplot(2,2,2);
stem(n, imag(x_n));  % Sanal k�s�m i�in 
title('Sanal K�s�m');
xlabel('n');
ylabel('Im{x[n]}');

subplot(2,2,3);
stem(n, abs(x_n));  % Genlik i�in 
title('Genlik');
xlabel('n');
ylabel('|x[n]|');

subplot(2,2,4);
stem(n, angle(x_n));  % Faz i�in 
title('Faz');
xlabel('n');
ylabel('?x[n]');



% 1 soru b ��kk�


figure;  
t = 0:0.01:250;  % S�rekli zaman aral���
x_t = exp(-(0.05 + 0.3j) * t);  % s�rekli zaman i�aretinin  denklemi


% b ��kk� i�in grafikler

subplot(2,2,1);
plot(t, real(x_t));
title('Ger�ek K�s�m');

subplot(2,2,2);
plot(t, imag(x_t));
title('Sanal K�s�m');

subplot(2,2,3);
plot(t, abs(x_t));
title('Genlik');

subplot(2,2,4);
plot(t, angle(x_t));
title('Faz');


%  2.soru

t = -8:0.0001:8;  % bir aral�k tan�mlama
x_t = zeros(size(t));  % x(t) i�aret vekt�r� 

% x(t) 
x_t(t>-1&t<0) = t(t>-1&t<0)+1;  % -1 < t < 0 aral���nda x(t) = t

x_t(t>=0&t<2) = 1;  % 0 <= t < 2 aral���nda x(t) = 1


% Orjinal x(t) ��areti:
% � -1<t<0 aral���nda x(t) = t + 1.
% � 0<t< 2 aral���nda x(t) = 1.
% � Di�er t�m t de�erleri i�in, x(t) = 0



% x(t) i�aretinin grafi�i
figure;
plot(t, x_t, 'LineWidth', 2);
title('x(t) ��areti');
xlabel('t');
ylabel('x(t)');



% 2.soru a ��kk�
x_2t = zeros(size(t));  % x(2t) i�in vekt�r 

% x(2t) 
x_2t(t > -0.5 & t < 0) = 2 * t(t > -0.5 & t < 0) + 1;  % -0.5 < t < 0 aral���nda x(2t) = 2t + 1
x_2t(t >= 0 & t < 1) = 1;  % 0 <= t < 1 aral���nda x(2t) = 1
% t < -0.5 ve t >= 1 i�in x(2t) zaten 0 olarak tan�ml�


% (2t) i�in aral�klar ��yledir: 
% -0.5 < t < 0 i�in, x(2t) = 2t+1.
% 0<t<1 i�in, a(2t) = 1. 
% t <-0.5 ve t > 1 i�in, x(2t) = 0.

% x(2t) i�aretini �izdir
figure;
plot(t, x_2t, 'LineWidth', 2);
title('x(2t) ��areti');
xlabel('t');
ylabel('x(2t)');
 




% 2.soru b ��kk� 

x_t_delayed = zeros(size(t));  % x(t - 5/2) i�in bir vekt�r

% x(t - 5/2) 
x_t_delayed(t > (-1 + 5/2) & t < (5/2)) = t(t > (-1 + 5/2) & t < (5/2)) - 5/2 + 1;  % -1 + 5/2 < t < 5/2 aral���nda x(t - 5/2) = t - 5/2 + 1
x_t_delayed(t >= (5/2) & t < (2 + 5/2)) = 1;  % 5/2 <= t < 2 + 5/2 aral���nda x(t - 5/2) = 1
% t < -1 + 5/2 ve t >= 2 + 5/2 i�in x(t - 5/2) zaten 0 olarak tan�ml�

% x(t - 5/2) i�aretinin grafi�i 
figure;
plot(t, x_t_delayed, 'LineWidth', 2);
title('x(t - 5/2) ��areti');
xlabel('t');
ylabel('x(t - 5/2)');


% x(t) i�aretini birim sa�a kayd�r�rsak, her aral�k birim sa�a kayacakt�r
% -1<t<0 aral��� birim sa�a kayd�r�ld���nda, -3.5< t-2.5 < -2.5 
% x(t)=t-2.5 + 1 = t - 1.5 olur.
% 0<t<2 aral��� birim sa�a kayd�r�ld���nda, -2.5 ? t <0.5  Bu aral�kta x(t-5/2)=1 olur.
% Dolay�s�yla, kayd�r�lm�� x(t-5/2) i�areti i�in aral�klar �unlard�r:
% -3.5<t<-2.5 i�in, x(t-5/2)=t-1.5
% -2.5<t<-0.5 i�in,  x(t-5/2)=1 
%  Di�er t�m t de�erleri i�in, x (t -5/2 ) = 0 d�r .




% 2.soru c ��kk� 

% x(-2t+5) 
x_neg2t_plus5 = zeros(size(t));  % x(-2t+5) i�aretinin vekt�r� 

% x(-2t+5) i�areti

t_prime = -2.*t+5;
x_neg2t_plus5(t_prime > -1 & t_prime < 0) = t_prime(t_prime > -1 & t_prime < 0) + 1;
x_neg2t_plus5(t_prime >= 0 & t_prime < 2) = 1;

% x(-2t+5) i�aretinin garfi�i 
figure;
plot(t, x_neg2t_plus5, 'LineWidth', 2);
title('x(-2t+5) ��areti');
xlabel('t');
ylabel('x(-2t+5)');


 
% tek ve �ift k�s�mlar�n�n� x(t) ve x(2t) i�in ayn� ��kmas�n� �l�ek olarak
% farkl� fakat g�r�n�� olarak ayn� geldigi g�zlemlendi. Beklenen bir durum
% g�zlemlenmi� oldu sadece �l�ekte bir oynama grafikte g�r�n�� olarak bi�e
% degi�tirmez sadece �l�ekte bir degi�iklik meydana gelir .
t = -10:0.1:10;


% x(t)  
x_1t = (t >= -1 & t < 0) .* (t + 1) + (t >= 0 & t < 2);   

% x(2t) 
x_2t = (2*t >= -1 & 2*t < 0) .* (2*t + 1) + (2*t >= 0 & 2*t < 2);  

% x(-2t) 
x_neg_2t = (-2*t >= -1 & -2*t < 0) .* (-2*t + 1) + (-2*t >= 0 & -2*t < 2);   

% x(t - 5/2) 
x_t_minus_52 = (t - 5/2 >= -1 & t - 5/2 < 0) .* (t - 5/2 + 1) + (t - 5/2 >= 0 & t - 5/2 < 2);  

% x(-t - 5/2) 
x_neg_t_minus_52 = (-t - 5/2 >= -1 & -t - 5/2 < 0) .* (-t - 5/2 + 1) + (-t - 5/2 >= 0 & -t - 5/2 < 2);  

% x(-2t + 5) 
x_neg_2t_plus_5 = (-2*t + 5 >= -1 & -2*t + 5 < 0) .* (-2*t + 5 + 1) + (-2*t + 5 >= 0 & -2*t + 5 < 2);  

% x(2t + 5) 
x_2t_plus_5 = (2*t + 5 >= -1 & 2*t + 5 < 0) .* (2*t + 5 + 1) + (2*t + 5 >= 0 & 2*t + 5 < 2);  

% odd and even 
x_2t_odd = (x_2t - x_neg_2t) * 0.5;  % x(2t) odd
x_2t_even = (x_2t + x_neg_2t) * 0.5;  % x(2t) even

x_t_minus_52_odd = (x_t_minus_52 - x_neg_t_minus_52) * 0.5;  % x(t - 5/2) odd
x_t_minus_52_even = (x_t_minus_52 + x_neg_t_minus_52) * 0.5;  % x(t - 5/2) even

x_neg_2t_plus_5_odd = (x_neg_2t_plus_5 - x_2t_plus_5) * 0.5;  % x(-2t + 5) odd
x_neg_2t_plus_5_even = (x_neg_2t_plus_5 + x_2t_plus_5) * 0.5;  % x(-2t + 5) even

% garfikeler


% x(2t) odd  Grafi�i
figure;
subplot(3,2,1);
plot(t, x_2t_odd);  
title('x(2t) Odd Component');
xlabel('t');
ylabel('Amplitude');


% x(2t) even  Grafi�i

subplot(3,2,2);
plot(t, x_2t_even);  
title('x(2t) Even Component');
xlabel('t');
ylabel('Amplitude');

% x(t - 5/2) odd grafi�i
subplot(3,2,3);
plot(t, x_t_minus_52_odd);  
title('x(t - 5/2) Odd Component');
xlabel('t');
ylabel('Amplitude');


% x(t - 5/2) even grafi�i
subplot(3,2,4);
plot(t, x_t_minus_52_even);  
title('x(t - 5/2) Even Component');
xlabel('t');
ylabel('Amplitude');


% x(-2t + 5) odd grafi�i
subplot(3,2,5);
plot(t, x_neg_2t_plus_5_odd);  
title('x(-2t + 5) Odd Component');
xlabel('t');
ylabel('Amplitude');


% x(-2t + 5) odd grafi�i
subplot(3,2,6);
plot(t, x_neg_2t_plus_5_even);  
title('x(-2t + 5) Even Component');
xlabel('t');
ylabel('Amplitude');










% aral�k belirleme 
n = -10:10 ;

%  x[n] and h[n]
x_n = ((1/4).^n) .* ((n >= 0) & (n <= 5));
h_n = ((1/2).^n) .* ((n >= -3) & (n <= 2));


y_n = zeros(1, length(n) + length(n) - 1);

%convolution 
for i = 1:length(n)
    for j = 1:length(n)
        if (i + j - 1 <= length(y_n))
            y_n(i + j - 1) = y_n(i + j - 1) + x_n(i) * h_n(j);
        end
    end
end

n_conv = (2 * min(n)):1:(2 * max(n));

% sinyalleri grafikleri

figure; 

stem(n, x_n);
title('x[n]');
xlabel('n');
ylabel('x[n]');
figure;

stem(n, h_n);
title('h[n]');
xlabel('n');
ylabel('h[n]');
figure;

stem(n_conv, y_n);
title('Convolution y[n] = x[n] * h[n]');
xlabel('n');
ylabel('y[n]');
t = -10:0.01:10;






% x(t) i�aretini �rnekleme (sampling)
% x(t) = sin(pi*t) i�in 0 <= t <= 2

x_t = sin(pi * t) .* (t >= 0 & t <= 2);


% h(t) = 1/2 i�in -2 <= t <= 5
h_t = 0.5 * (t >= -2 & t <= 5);
y_t = zeros(1, length(t) + length(t) - 1);



% Konvol�syon 
for i = 1:length(t)
    for j = 1:length(t)
        if (i + j - 1 <= length(y_t))
            y_t(i + j - 1) = y_t(i + j - 1) + x_t(i) * h_t(j) * 0.01;  % Delta t �arpan�n� eklemeyi unutmay�n
        end
    end
end


% x(t) Grafi�i
figure;
plot(t, x_t);
title('x(t) = sin(pi * t) * (u(t) - u(t - 2))');
xlabel('t');
ylabel('x(t)');
grid on ;
% h(t)  Grafi�i
figure;
plot(t, h_t);
title('h(t) = 0.5 * (u(t + 2) - u(t - 5))');
xlabel('t');
ylabel('h(t)');
grid on ;
% Convolution y(t) = x(t) * h(t)
figure;
plot(-20:0.01:20, y_t);
title('Convolution y(t) = x(t) * h(t)');
xlabel('t');
ylabel('y(t)');
grid on ;




