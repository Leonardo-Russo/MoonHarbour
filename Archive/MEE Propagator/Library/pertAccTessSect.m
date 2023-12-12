function [ar, aE, aN] = pertAccTessSect(rr, lat, lamg)


global Rm muM
global J2_2 J3_1 J3_2 J3_3 J4_1 J4_2 J5_1 J6_1 J7_1 J9_1 J12_1;
global lamg22 lamg31 lamg32 lamg33 lamg41 lamg42 lamg51 lamg61 lamg71 lamg91 lamg121;

arJ2_2=        -(9*J2_2*Rm^2*muM*cos(2*lamg - 2*lamg22)*cos(lat)^2)/rr^4;
aEJ2_2=          -(6*J2_2*Rm^2*muM*sin(2*lamg - 2*lamg22)*cos(lat))/rr^4;
aNJ2_2= -(6*J2_2*Rm^2*muM*cos(2*lamg - 2*lamg22)*cos(lat)*sin(lat))/rr^4;

arJ3_1=                       -(4*J3_1*Rm^3*muM*cos(lamg - lamg31)*(cos(lat) - 5*cos(lat)*((3*cos(lat)^2)/2 - 1)))/rr^5;
aEJ3_1=              -(J3_1*Rm^3*muM*sin(lamg - lamg31)*(cos(lat) - 5*cos(lat)*((3*cos(lat)^2)/2 - 1)))/(rr^5*cos(lat));
aNJ3_1= (J3_1*Rm^3*muM*cos(lamg - lamg31)*(5*sin(lat)*((3*cos(lat)^2)/2 - 1) - sin(lat) + 15*cos(lat)^2*sin(lat)))/rr^5;

% arJ3_2=                                                               -(60*J3_2*R^3*mi*cos(2*lamg - 2*lamg32)*cos(lat)^2*sin(lat))/rr^5;
% aEJ3_2=                                                                 -(30*J3_2*R^3*mi*sin(2*lamg - 2*lamg32)*cos(lat)*sin(lat))/rr^5;
% aNJ3_2= ((15*J3_2*R^3*mi*cos(2*lamg - 2*lamg32)*cos(lat)^3)/rr^4 - (30*J3_2*R^3*mi*cos(2*lamg - 2*lamg32)*cos(lat)*sin(lat)^2)/rr^4)/rr;
% 
% arJ3_3=           -(60*J3_3*R^3*mi*cos(3*lamg - 3*lamg33)*cos(lat)^3)/rr^5;
% aEJ3_3=          -(45*J3_3*R^3*mi*sin(3*lamg - 3*lamg33)*cos(lat)^2)/rr^5;
% aNJ3_3= -(45*J3_3*R^3*mi*cos(3*lamg - 3*lamg33)*cos(lat)^2*sin(lat))/rr^5;
%  
% arJ4_1=                                                                -(5*J4_1*R^4*mi*cos(lamg - lamg41)*(3*cos(lat)*sin(lat) + (7*cos(lat)*sin(lat)*(5*sin(lat)^2 - 3))/2))/rr^6;
% aEJ4_1=                                                       -(J4_1*R^4*mi*sin(lamg - lamg41)*(3*cos(lat)*sin(lat) + (7*cos(lat)*sin(lat)*(5*sin(lat)^2 - 3))/2))/(rr^6*cos(lat));
% aNJ4_1= (J4_1*R^4*mi*cos(lamg - lamg41)*(3*cos(lat)^2 - 3*sin(lat)^2 + (7*cos(lat)^2*(5*sin(lat)^2 - 3))/2 - (7*sin(lat)^2*(5*sin(lat)^2 - 3))/2 + 35*cos(lat)^2*sin(lat)^2))/rr^6;
% 
% arJ4_2=                                                                                           -(5*J4_2*R^4*mi*cos(2*lamg - 2*lamg42)*(3*cos(lat)^2 + 7*cos(lat)*(cos(lat) - 5*cos(lat)*((3*cos(lat)^2)/2 - 1))))/rr^6;
% aEJ4_2=                                                                                -(2*J4_2*R^4*mi*sin(2*lamg - 2*lamg42)*(3*cos(lat)^2 + 7*cos(lat)*(cos(lat) - 5*cos(lat)*((3*cos(lat)^2)/2 - 1))))/(rr^6*cos(lat));
% aNJ4_2= -(J4_2*R^4*mi*cos(2*lamg - 2*lamg42)*(6*cos(lat)*sin(lat) - 7*cos(lat)*(5*sin(lat)*((3*cos(lat)^2)/2 - 1) - sin(lat) + 15*cos(lat)^2*sin(lat)) + 7*sin(lat)*(cos(lat) - 5*cos(lat)*((3*cos(lat)^2)/2 - 1))))/rr^6;
% 
% arJ5_1=                                                                                         -(6*J5_1*R^5*mi*cos(lamg - lamg51)*(cos(lat) - 5*cos(lat)*((3*cos(lat)^2)/2 - 1) + 9*cos(lat)*((35*cos(lat)^4)/8 - 5*cos(lat)^2 + 1)))/rr^7;
% aEJ5_1=                                                                               -(J5_1*R^5*mi*sin(lamg - lamg51)*(cos(lat) - 5*cos(lat)*((3*cos(lat)^2)/2 - 1) + 9*cos(lat)*((35*cos(lat)^4)/8 - 5*cos(lat)^2 + 1)))/(rr^7*cos(lat));
% aNJ5_1= (J5_1*R^5*mi*cos(lamg - lamg51)*(9*cos(lat)*(10*cos(lat)*sin(lat) - (35*cos(lat)^3*sin(lat))/2) - sin(lat) + 5*sin(lat)*((3*cos(lat)^2)/2 - 1) + 15*cos(lat)^2*sin(lat) - 9*sin(lat)*((35*cos(lat)^4)/8 - 5*cos(lat)^2 + 1)))/rr^7;
%  
% arJ6_1=                                                                                                                                                                                                -(7*J6_1*R^6*mi*cos(lamg - lamg61)*(3*cos(lat)*sin(lat) + (7*cos(lat)*sin(lat)*(5*sin(lat)^2 - 3))/2 + (11*cos(lat)*sin(lat)*(63*sin(lat)^4 - 70*sin(lat)^2 + 15))/8))/rr^8;
% aEJ6_1=                                                                                                                                                                                      -(J6_1*R^6*mi*sin(lamg - lamg61)*(3*cos(lat)*sin(lat) + (7*cos(lat)*sin(lat)*(5*sin(lat)^2 - 3))/2 + (11*cos(lat)*sin(lat)*(63*sin(lat)^4 - 70*sin(lat)^2 + 15))/8))/(rr^8*cos(lat));
% aNJ6_1= (J6_1*R^6*mi*cos(lamg - lamg61)*((11*cos(lat)^2*(63*sin(lat)^4 - 70*sin(lat)^2 + 15))/8 + 3*cos(lat)^2 - (11*sin(lat)^2*(63*sin(lat)^4 - 70*sin(lat)^2 + 15))/8 - 3*sin(lat)^2 + (7*cos(lat)^2*(5*sin(lat)^2 - 3))/2 - (7*sin(lat)^2*(5*sin(lat)^2 - 3))/2 + 35*cos(lat)^2*sin(lat)^2 - (11*cos(lat)*sin(lat)*(140*cos(lat)*sin(lat) - 252*cos(lat)*sin(lat)^3))/8))/rr^8;
%  
% arJ7_1=                                                                                                                                                                                                                                                                                                                                                  (7*J7_1*R^7*mi*cos(lat)*cos(lamg - lamg71)*(429*cos(lat)^6 - 792*cos(lat)^4 + 432*cos(lat)^2 - 64))/(2*rr^9);
% aEJ7_1=                                                                                                                                                                                                                                                                                                                                                         (7*J7_1*R^7*mi*sin(lamg - lamg71)*(429*cos(lat)^6 - 792*cos(lat)^4 + 432*cos(lat)^2 - 64))/(16*rr^9);
% aNJ7_1= (J7_1*R^7*mi*(2*sin(lamg/2 - lamg71/2)^2 - 1)*(sin(lat)*(- (819*sin(2*lat)^2)/128 + (3003*sin(3*lat)^2)/256 + (1365*sin(lat)^2)/256 - 65/16) - (11*sin(lat))/4 - (15*sin(3*lat))/4 + sin(lat)*(- (315*sin(2*lat)^2)/32 + (45*sin(lat)^2)/8 + 27/8) + sin(lat)*((15*sin(lat)^2)/2 - 5/2) + (2*sin(lat/2)^2 - 1)*((45*sin(2*lat))/8 - (315*sin(4*lat))/16) + (2*sin(lat/2)^2 - 1)*((1365*sin(2*lat))/256 - (819*sin(4*lat))/64 + (9009*sin(6*lat))/256)))/rr^9;
%  
% arJ9_1=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          -(225*J9_1*R^9*mi*cos(lat)*cos(lamg - lamg91)*(2431*cos(lat)^8 - 5720*cos(lat)^6 + 4576*cos(lat)^4 - 1408*cos(lat)^2 + 128))/(64*rr^11);
% aEJ9_1=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   -(45*J9_1*R^9*mi*sin(lamg - lamg91)*(2431*cos(lat)^8 - 5720*cos(lat)^6 + 4576*cos(lat)^4 - 1408*cos(lat)^2 + 128))/(128*rr^11);
% aNJ9_1= (J9_1*R^9*mi*(2*sin(lamg/2 - lamg91/2)^2 - 1)*(sin(lat)*(- (11781*sin(2*lat)^2)/2048 + (7293*sin(3*lat)^2)/1024 - (109395*sin(4*lat)^2)/8192 + (5355*sin(lat)^2)/1024 + 595/128) - (11*sin(lat))/4 - (15*sin(3*lat))/4 + sin(lat)*(- (819*sin(2*lat)^2)/128 + (3003*sin(3*lat)^2)/256 + (1365*sin(lat)^2)/256 - 65/16) + sin(lat)*(- (315*sin(2*lat)^2)/32 + (45*sin(lat)^2)/8 + 27/8) + sin(lat)*((15*sin(lat)^2)/2 - 5/2) + (2*sin(lat/2)^2 - 1)*((45*sin(2*lat))/8 - (315*sin(4*lat))/16) + (2*sin(lat/2)^2 - 1)*((5355*sin(2*lat))/1024 - (11781*sin(4*lat))/1024 + (21879*sin(6*lat))/1024 - (109395*sin(8*lat))/2048) + (2*sin(lat/2)^2 - 1)*((1365*sin(2*lat))/256 - (819*sin(4*lat))/64 + (9009*sin(6*lat))/256)))/rr^11;
%  
% arJ12_1=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     -(13*J12_1*R^12*mi*cos(lamg - lamg121)*(3*cos(lat)*sin(lat) + (7*cos(lat)*sin(lat)*(5*sin(lat)^2 - 3))/2 + (15*cos(lat)*sin(lat)*(429*sin(lat)^6 - 693*sin(lat)^4 + 315*sin(lat)^2 - 35))/16 + (23*cos(lat)*sin(lat)*(88179*sin(lat)^10 - 230945*sin(lat)^8 + 218790*sin(lat)^6 - 90090*sin(lat)^4 + 15015*sin(lat)^2 - 693))/256 + (11*cos(lat)*sin(lat)*(63*sin(lat)^4 - 70*sin(lat)^2 + 15))/8 + (19*cos(lat)*sin(lat)*(12155*sin(lat)^8 - 25740*sin(lat)^6 + 18018*sin(lat)^4 - 4620*sin(lat)^2 + 315))/128))/rr^14;
% aEJ12_1=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             -(J12_1*R^12*mi*sin(lamg - lamg121)*(3*sin(lat) + (19*sin(lat)*(12155*sin(lat)^8 - 25740*sin(lat)^6 + 18018*sin(lat)^4 - 4620*sin(lat)^2 + 315))/128 + (7*sin(lat)*(5*sin(lat)^2 - 3))/2 + (15*sin(lat)*(429*sin(lat)^6 - 693*sin(lat)^4 + 315*sin(lat)^2 - 35))/16 + (23*sin(lat)*(88179*sin(lat)^10 - 230945*sin(lat)^8 + 218790*sin(lat)^6 - 90090*sin(lat)^4 + 15015*sin(lat)^2 - 693))/256 + (11*sin(lat)*(63*sin(lat)^4 - 70*sin(lat)^2 + 15))/8))/rr^14;
% aNJ12_1= (J12_1*R^12*mi*cos(lamg - lamg121)*(3*cos(2*lat) - (35*cos(4*lat))/8 + cos(2*lat)*((77*cos(2*lat))/16 + (693*cos(4*lat))/64 + 319/64) + (sin(2*lat)*((877565*sin(2*lat))/32768 + (133055*sin(4*lat))/4096 + (5413395*sin(6*lat))/65536 + (482885*sin(8*lat))/8192 + (10140585*sin(10*lat))/65536))/2 + cos(2*lat)*((14839*cos(2*lat))/2048 + (51623*cos(4*lat))/4096 + (13585*cos(6*lat))/2048 + (230945*cos(8*lat))/16384 + 101251/16384) - (sin(2*lat)*((77*sin(2*lat))/8 + (693*sin(4*lat))/16))/2 - (sin(2*lat)*((14839*sin(2*lat))/1024 + (51623*sin(4*lat))/1024 + (40755*sin(6*lat))/1024 + (230945*sin(8*lat))/2048))/2 - cos(2*lat)*((35*cos(2*lat))/4 + 7/4) - cos(2*lat)*((5805*cos(2*lat))/512 + (1485*cos(4*lat))/256 + (6435*cos(6*lat))/512 + 795/256) - cos(2*lat)*((877565*cos(2*lat))/65536 + (133055*cos(4*lat))/16384 + (1804465*cos(6*lat))/131072 + (482885*cos(8*lat))/65536 + (2028117*cos(10*lat))/131072 + 271423/65536) + (sin(2*lat)*((5805*sin(2*lat))/256 + (1485*sin(4*lat))/64 + (19305*sin(6*lat))/256))/2 + 35/8))/rr^14;
 


%%% expressions evaluated through vpa(...,22)

% arJ3_1=                           -(4.0*J3_1*R^3*mi*cos(lamg - 1.0*lamg31)*(cos(lat) - 5.0*cos(lat)*(1.5*cos(lat)^2 - 1.0)))/rr^5;
% aEJ3_1=                -(1.0*J3_1*R^3*mi*sin(lamg - 1.0*lamg31)*(cos(lat) - 5.0*cos(lat)*(1.5*cos(lat)^2 - 1.0)))/(rr^5*cos(lat));
% aNJ3_1= (J3_1*R^3*mi*cos(lamg - 1.0*lamg31)*(5.0*sin(lat)*(1.5*cos(lat)^2 - 1.0) - 1.0*sin(lat) + 15.0*cos(lat)^2*sin(lat)))/rr^5;

arJ3_2=                                                                     -(60.0*J3_2*Rm^3*muM*cos(2.0*lamg - 2.0*lamg32)*cos(lat)^2*sin(lat))/rr^5;
aEJ3_2=                                                                       -(30.0*J3_2*Rm^3*muM*sin(2.0*lamg - 2.0*lamg32)*cos(lat)*sin(lat))/rr^5;
aNJ3_2= ((15.0*J3_2*Rm^3*muM*cos(2.0*lamg - 2.0*lamg32)*cos(lat)^3)/rr^4 - (30.0*J3_2*Rm^3*muM*cos(2.0*lamg - 2.0*lamg32)*cos(lat)*sin(lat)^2)/rr^4)/rr;
                                                                                                                                                                                                                                                                         
arJ3_3=          -(60.0*J3_3*Rm^3*muM*cos(3.0*lamg - 3.0*lamg33)*cos(lat)^3)/rr^5;
aEJ3_3=          -(45.0*J3_3*Rm^3*muM*sin(3.0*lamg - 3.0*lamg33)*cos(lat)^2)/rr^5;
aNJ3_3=-(45.0*J3_3*Rm^3*muM*cos(3.0*lamg - 3.0*lamg33)*cos(lat)^2*sin(lat))/rr^5;
                                                                                                                                                                                                            
arJ4_1=                                                                    -(5.0*J4_1*Rm^4*muM*cos(lamg - 1.0*lamg41)*(3.0*cos(lat)*sin(lat) + 3.5*cos(lat)*sin(lat)*(5.0*sin(lat)^2 - 3.0)))/rr^6;
aEJ4_1=                                                         -(1.0*J4_1*Rm^4*muM*sin(lamg - 1.0*lamg41)*(3.0*cos(lat)*sin(lat) + 3.5*cos(lat)*sin(lat)*(5.0*sin(lat)^2 - 3.0)))/(rr^6*cos(lat));
aNJ4_1=(J4_1*Rm^4*muM*cos(lamg - 1.0*lamg41)*(3.0*cos(lat)^2 - 3.0*sin(lat)^2 + 3.5*cos(lat)^2*(5.0*sin(lat)^2 - 3.0) - 3.5*sin(lat)^2*(5.0*sin(lat)^2 - 3.0) + 35.0*cos(lat)^2*sin(lat)^2))/rr^6;
                                                                                          
arJ4_2=                                                                                                       -(5.0*J4_2*Rm^4*muM*cos(2.0*lamg - 2.0*lamg42)*(3.0*cos(lat)^2 + 7.0*cos(lat)*(cos(lat) - 5.0*cos(lat)*(1.5*cos(lat)^2 - 1.0))))/rr^6;
aEJ4_2=                                                                                            -(2.0*J4_2*Rm^4*muM*sin(2.0*lamg - 2.0*lamg42)*(3.0*cos(lat)^2 + 7.0*cos(lat)*(cos(lat) - 5.0*cos(lat)*(1.5*cos(lat)^2 - 1.0))))/(rr^6*cos(lat));
aNJ4_2= -(1.0*J4_2*Rm^4*muM*cos(2.0*lamg - 2.0*lamg42)*(6.0*cos(lat)*sin(lat) - 7.0*cos(lat)*(5.0*sin(lat)*(1.5*cos(lat)^2 - 1.0) - 1.0*sin(lat) + 15.0*cos(lat)^2*sin(lat)) + 7.0*sin(lat)*(cos(lat) - 5.0*cos(lat)*(1.5*cos(lat)^2 - 1.0))))/rr^6;
                                                   
arJ5_1=                                                                                              -(6.0*J5_1*Rm^5*muM*cos(lamg - 1.0*lamg51)*(cos(lat) - 5.0*cos(lat)*(1.5*cos(lat)^2 - 1.0) + 9.0*cos(lat)*(4.375*cos(lat)^4 - 5.0*cos(lat)^2 + 1.0)))/rr^7;
aEJ5_1=                                                                                   -(1.0*J5_1*Rm^5*muM*sin(lamg - 1.0*lamg51)*(cos(lat) - 5.0*cos(lat)*(1.5*cos(lat)^2 - 1.0) + 9.0*cos(lat)*(4.375*cos(lat)^4 - 5.0*cos(lat)^2 + 1.0)))/(rr^7*cos(lat));
aNJ5_1=(J5_1*Rm^5*muM*cos(lamg - 1.0*lamg51)*(9.0*cos(lat)*(10.0*cos(lat)*sin(lat) - 17.5*cos(lat)^3*sin(lat)) - 1.0*sin(lat) + 5.0*sin(lat)*(1.5*cos(lat)^2 - 1.0) + 15.0*cos(lat)^2*sin(lat) - 9.0*sin(lat)*(4.375*cos(lat)^4 - 5.0*cos(lat)^2 + 1.0)))/rr^7;
             
arJ6_1=                                                                                                                                                                                                           -(7.0*J6_1*Rm^6*muM*cos(lamg - 1.0*lamg61)*(3.0*cos(lat)*sin(lat) + 1.375*cos(lat)*sin(lat)*(63.0*sin(lat)^4 - 70.0*sin(lat)^2 + 15.0) + 3.5*cos(lat)*sin(lat)*(5.0*sin(lat)^2 - 3.0)))/rr^8;
aEJ6_1=                                                                                                                                                                                                -(1.0*J6_1*Rm^6*muM*sin(lamg - 1.0*lamg61)*(3.0*cos(lat)*sin(lat) + 1.375*cos(lat)*sin(lat)*(63.0*sin(lat)^4 - 70.0*sin(lat)^2 + 15.0) + 3.5*cos(lat)*sin(lat)*(5.0*sin(lat)^2 - 3.0)))/(rr^8*cos(lat));
aNJ6_1= (J6_1*Rm^6*muM*cos(lamg - 1.0*lamg61)*(3.0*cos(lat)^2 - 3.0*sin(lat)^2 + 3.5*cos(lat)^2*(5.0*sin(lat)^2 - 3.0) - 3.5*sin(lat)^2*(5.0*sin(lat)^2 - 3.0) + 35.0*cos(lat)^2*sin(lat)^2 + 1.375*cos(lat)^2*(63.0*sin(lat)^4 - 70.0*sin(lat)^2 + 15.0) - 1.375*sin(lat)^2*(63.0*sin(lat)^4 - 70.0*sin(lat)^2 + 15.0) - 1.375*cos(lat)*sin(lat)*(140.0*cos(lat)*sin(lat) - 252.0*cos(lat)*sin(lat)^3)))/rr^8;

arJ7_1=       (3.5*J7_1*Rm^7*muM*cos(lamg - 1.0*lamg71)*cos(lat)*(432.0*cos(lat)^2 - 792.0*cos(lat)^4 + 429.0*cos(lat)^6 - 64.0))/rr^9;
aEJ7_1=             -(0.4375*J7_1*Rm^7*muM*sin(lamg - 1.0*lamg71)*(135.0*sin(lat)^2 - 495.0*sin(lat)^4 + 429.0*sin(lat)^6 - 5.0))/rr^9;
aNJ7_1= (0.4375*J7_1*Rm^7*muM*cos(lamg - 1.0*lamg71)*sin(lat)*(1296.0*cos(lat)^2 - 3960.0*cos(lat)^4 + 3003.0*cos(lat)^6 - 64.0))/rr^9;

arJ9_1=     -(3.515625*J9_1*Rm^9*muM*cos(lamg - 1.0*lamg91)*cos(lat)*(4576.0*cos(lat)^4 - 1408.0*cos(lat)^2 - 5720.0*cos(lat)^6 + 2431.0*cos(lat)^8 + 128.0))/rr^11;
aEJ9_1=                -(0.3515625*J9_1*Rm^9*muM*sin(lamg - 1.0*lamg91)*(2002.0*sin(lat)^4 - 308.0*sin(lat)^2 - 4004.0*sin(lat)^6 + 2431.0*sin(lat)^8 + 7.0))/rr^11;
aNJ9_1=-(0.3515625*J9_1*Rm^9*muM*cos(lamg - 1.0*lamg91)*sin(lat)*(22880.0*cos(lat)^4 - 4224.0*cos(lat)^2 - 40040.0*cos(lat)^6 + 21879.0*cos(lat)^8 + 128.0))/rr^11;

arJ12_1=             (1.98046875*J12_1*Rm^12*muM*cos(lamg - 1.0*lamg121)*cos(lat)*sin(lat)*(9600.0*cos(lat)^2 - 54400.0*cos(lat)^4 + 129200.0*cos(lat)^6 - 135660.0*cos(lat)^8 + 52003.0*cos(lat)^10 - 512.0))/rr^14;
aEJ12_1=                     -(0.15234375*J12_1*Rm^12*muM*sin(lamg - 1.0*lamg121)*sin(lat)*(5775.0*sin(lat)^2 - 39270.0*sin(lat)^4 + 106590.0*sin(lat)^6 - 124355.0*sin(lat)^8 + 52003.0*sin(lat)^10 - 231.0))/rr^14;
aNJ12_1= -(0.15234375*J12_1*Rm^12*muM*cos(lamg - 1.0*lamg121)*(310400.0*cos(lat)^4 - 29824.0*cos(lat)^2 - 1230800.0*cos(lat)^6 + 2254540.0*cos(lat)^8 - 1928633.0*cos(lat)^10 + 624036.0*cos(lat)^12 + 512.0))/rr^14;

ar=arJ2_2+arJ3_1+arJ3_2+arJ3_3+arJ4_1+arJ4_2+arJ5_1+arJ6_1+arJ7_1+arJ9_1+arJ12_1;
aE=aEJ2_2+aEJ3_1+aEJ3_2+aEJ3_3+aEJ4_1+aEJ4_2+aEJ5_1+aEJ6_1+aEJ7_1+aEJ9_1+aEJ12_1;
aN=aNJ2_2+aNJ3_1+aNJ3_2+aNJ3_3+aNJ4_1+aNJ4_2+aNJ5_1+aNJ6_1+aNJ7_1+aNJ9_1+aNJ12_1;

%[aNJ2_2,aNJ3_1,aNJ3_2,aNJ3_3,aNJ4_1,aNJ4_2,aNJ5_1,aNJ6_1,aNJ7_1,aNJ9_1,aNJ12_1]
%[arJ2,arJ3,arJ9,arJ10,arJ11,arJ12,arJ17,arJ22,arJ23,arJ24,arJ25,arJ27,arJ29,arJ39,arJ46,arJ54,arJ67,arJ69]


end