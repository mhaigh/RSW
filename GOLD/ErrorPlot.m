% plot the error
clear all
close all

error0=[1.532644579525016e-12,4.893713456569950e-10,2.204910840777935e-05];
error1=[9.491555231854979e-14,1.717731194865110e-11,1.178135813865897e-06];
error2=[2.222199076963427e-14,3.526041801637414e-12,1.586021931459977e-07];
error3=[2.840098546457533e-15,4.372185150410833e-13,1.787922495176188e-08];
error4=[2.264277112121863e-16,3.448203860064709e-14,1.495237945964300e-09];
error5=[1.516989335038468e-17,2.299371259414864e-15,1.281823837212829e-10];
error6=[9.667619396422008e-19,1.462122593411069e-16,1.520317343447664e-11];
error7=[6.077933502675291e-20,9.183089835515520e-18,2.692904019037603e-12];

res=[2^5,2^6,2^7,2^8,2^9,2^10,2^11];
uerrors=[error1(1),error2(1),error3(1),error4(1),error5(1),error6(1),error7(1)];

figure(1);
ax1=subplot(1,1,1);
scatter(ax1,res,uerrors,'x');...
    title('u momentum equation error'); xlabel('Grid Resolution'); ylabel('Mean Square Error');...
    set(gca,'XTick',[2^5,2^6,2^7,2^8,2^9,2^10,2^11]); axis([2^5 2^11 10e-21 10e-13]);...
    set(gca,'xscale','log','yscale','log'); 




figure(2);
loglog([2^5,2^6,2^7,2^8,2^9,2^10,2^11],[error1(2),error2(2),error3(2),error4(2),error5(2),error6(2),error7(2)],'-x');...
    title('v momentum equation error'); xlabel('Grid Resolution'); ylabel('Mean Square Error');...
    set(gca,'XTick',[2^5,2^6,2^7,2^8,2^9,2^10,2^11]); axis([2^5 2^11 10e-18 10e-10]);

figure(3);
loglog([2^5,2^6,2^7,2^8,2^9,2^10,2^11],[error1(3),error2(3),error3(3),error4(3),error5(3),error6(3),error7(3)],'-x');...
    title('h continuity equation error'); xlabel('Grid Resolution'); ylabel('Mean Square Error');...
    set(gca,'XTick',[2^5,2^6,2^7,2^8,2^9,2^10,2^11]); axis([2^5 2^11 10e-13 10e-6]);