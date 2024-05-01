clear;
figure;
% ****** define the physical paramter ******
Lx = 64;
U = 8;
V = -1.25;
hole = 16; % the number of hole, = doping * (2 * Lx)
t2 = -0.3; % actually won't be used in file name

% ****** accuracy parameter ******
D = 500; % bond dimension
guass_factor = 1/300;

% ****** define the figure range ****** %
delta_kx = 0.01;
k_set = -pi:delta_kx:pi;
omega_min = 0;
omega_max = 1.5;
omega_set = omega_max:-0.02:omega_min;

omega_size = numel(omega_set);

%_dynamichubbardLx64U8V-1.25hole8D500.json
FileNamePostfix=['_dynamichubbardLx',num2str(Lx),'U',num2str(U),...
    'V',num2str(V),'hole',num2str(hole),'D',num2str(D),'.json'];


[time, x_set, G_t_x0] =  ReadSpinCorrData(Lx, 0, FileNamePostfix);
[~, ~, G_t_x1] =  ReadSpinCorrData(Lx, 1, FileNamePostfix);
[~, ~, G_t_x2] =  ReadSpinCorrData(Lx, 2, FileNamePostfix);

[A0,Api] =  CalA_k_omega( k_set, omega_set, time, x_set, G_t_x0, guass_factor);
   % + 0.5 * CalA_k_omega( k_set, omega_set, time, x_set, G_t_x1, guass_factor) ...
   % + 0.5 * CalA_k_omega( k_set, omega_set, time, x_set, G_t_x2, guass_factor);


colormap(jet);
colormap parula;
imagesc(k_set, omega_set, A0); hold on;
colorbar
set(gca,'YDir','normal');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
%set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$k$','Interpreter','latex');
ylabel('$\omega/t$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

set(gca,'linewidth',1.5);
set(gcf,'position',[1000,1000,450,400]);
set(gca, 'XTick', [0,pi/2,pi]);
set(gca,'XTickLabel',{'0','\pi/2','\pi',});


figure;


colormap(jet);
colormap parula;
imagesc(k_set, omega_set, Api); hold on;
colorbar
set(gca,'YDir','normal');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
%set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$k$','Interpreter','latex');
ylabel('$\omega/t$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

set(gca,'linewidth',1.5);
set(gcf,'position',[1000,1000,450,400]);
set(gca, 'XTick', [0,pi/2,pi]);
set(gca,'XTickLabel',{'0','\pi/2','\pi',});


% set(gca,'Colormap',...
%     [1 1 1;0.987653017663562 0.989735773691646 0.99943037532401;0.975306035327124 0.979471547383293 0.99886075064802;0.962959052990686 0.969207321074939 0.998291125972031;0.950612070654248 0.958943094766586 0.997721501296041;0.93826508831781 0.948678868458232 0.997151876620051;0.925918105981372 0.938414642149879 0.996582251944061;0.913571123644934 0.928150415841525 0.996012627268071;0.901224141308495 0.917886189533172 0.995443002592081;0.888877158972057 0.907621963224818 0.994873377916091;0.876530176635619 0.897357736916465 0.994303753240102;0.864183194299181 0.887093510608111 0.993734128564112;0.851836211962743 0.876829284299758 0.993164503888122;0.839489229626305 0.866565057991404 0.992594879212132;0.827142247289867 0.856300831683051 0.992025254536142;0.814795264953429 0.846036605374697 0.991455629860152;0.802448282616991 0.835772379066344 0.990886005184163;0.790101300280553 0.82550815275799 0.990316380508173;0.777754317944115 0.815243926449637 0.989746755832183;0.765407335607677 0.804979700141283 0.989177131156193;0.753060353271239 0.79471547383293 0.988607506480203;0.740713370934801 0.784451247524576 0.988037881804213;0.728366388598362 0.774187021216223 0.987468257128224;0.716019406261924 0.763922794907869 0.986898632452234;0.703672423925486 0.753658568599515 0.986329007776244;0.691325441589048 0.743394342291162 0.985759383100254;0.67897845925261 0.733130115982808 0.985189758424264;0.666631476916172 0.722865889674455 0.984620133748274;0.654284494579734 0.712601663366101 0.984050509072284;0.641937512243296 0.702337437057748 0.983480884396295;0.629590529906858 0.692073210749394 0.982911259720305;0.61724354757042 0.681808984441041 0.982341635044315;0.604896565233982 0.671544758132687 0.981772010368325;0.592549582897544 0.661280531824334 0.981202385692335;0.580202600561106 0.65101630551598 0.980632761016345;0.567855618224667 0.640752079207627 0.980063136340356;0.55550863588823 0.630487852899273 0.979493511664366;0.543161653551791 0.62022362659092 0.978923886988376;0.530814671215353 0.609959400282566 0.978354262312386;0.518467688878915 0.599695173974213 0.977784637636396;0.506120706542477 0.589430947665859 0.977215012960406;0.493773724206039 0.579166721357506 0.976645388284416;0.481426741869601 0.568902495049152 0.976075763608427;0.469079759533163 0.558638268740798 0.975506138932437;0.456732777196725 0.548374042432445 0.974936514256447;0.444385794860287 0.538109816124091 0.974366889580457;0.432038812523849 0.527845589815738 0.973797264904467;0.419691830187411 0.517581363507384 0.973227640228478;0.407344847850973 0.507317137199031 0.972658015552488;0.394997865514535 0.497052910890677 0.972088390876498;0.382650883178097 0.486788684582324 0.971518766200508;0.370303900841658 0.47652445827397 0.970949141524518;0.35795691850522 0.466260231965617 0.970379516848528;0.345609936168782 0.455996005657263 0.969809892172538;0.333262953832344 0.44573177934891 0.969240267496549;0.320915971495906 0.435467553040556 0.968670642820559;0.308568989159468 0.425203326732203 0.968101018144569;0.29622200682303 0.414939100423849 0.967531393468579;0.288685559184855 0.414653307421329 0.965548371097328;0.281149111546679 0.414367514418808 0.963565348726077;0.273612663908504 0.414081721416287 0.961582326354827;0.266076216270329 0.413795928413767 0.959599303983576;0.258539768632153 0.413510135411246 0.957616281612325;0.251003320993978 0.413224342408726 0.955633259241074;0.243466873355802 0.412938549406205 0.953650236869823;0.235930425717627 0.412652756403685 0.951667214498573;0.228393978079452 0.412366963401164 0.949684192127322;0.220857530441276 0.412081170398643 0.947701169756071;0.213321082803101 0.411795377396123 0.94571814738482;0.205784635164926 0.411509584393602 0.943735125013569;0.19824818752675 0.411223791391082 0.941752102642319;0.190711739888575 0.410937998388561 0.939769080271068;0.183175292250399 0.41065220538604 0.937786057899817;0.179308950812965 0.41797926880105 0.934724710566611;0.175442609375532 0.42530633221606 0.931663363233404;0.171576267938098 0.43263339563107 0.928602015900197;0.167709926500664 0.43996045904608 0.925540668566991;0.16384358506323 0.44728752246109 0.922479321233784;0.159977243625796 0.4546145858761 0.919417973900578;0.156110902188362 0.46194164929111 0.916356626567371;0.152244560750928 0.46926871270612 0.913295279234165;0.148378219313494 0.47659577612113 0.910233931900958;0.14451187787606 0.48392283953614 0.907172584567752;0.140645536438626 0.49124990295115 0.904111237234545;0.136779195001192 0.49857696636616 0.901049889901339;0.132912853563758 0.50590402978117 0.897988542568132;0.129046512126324 0.51323109319618 0.894927195234925;0.12518017068889 0.52055815661119 0.891865847901719;0.121313829251457 0.5278852200262 0.888804500568512;0.117447487814023 0.53521228344121 0.885743153235306;0.113581146376589 0.54253934685622 0.882681805902099;0.109714804939155 0.54986641027123 0.879620458568893;0.105848463501721 0.55719347368624 0.876559111235686;0.101982122064287 0.564520537101249 0.87349776390248;0.098115780626853 0.571847600516259 0.870436416569273;0.094249439189419 0.579174663931269 0.867375069236067;0.0903830977519851 0.586501727346279 0.86431372190286;0.0865167563145512 0.593828790761289 0.861252374569654;0.0826504148771173 0.601155854176299 0.858191027236447;0.0787840734396833 0.608482917591309 0.85512967990324;0.0749177320022494 0.615809981006319 0.852068332570034;0.0710513905648155 0.623137044421329 0.849006985236827;0.0671850491273815 0.630464107836339 0.845945637903621;0.0633187076899476 0.637791171251349 0.842884290570414;0.0594523662525137 0.645118234666359 0.839822943237208;0.0555860248150798 0.652445298081369 0.836761595904001;0.0517196833776458 0.659772361496379 0.833700248570795;0.0478533419402119 0.667099424911389 0.830638901237588;0.043987000502778 0.674426488326399 0.827577553904382;0.040120659065344 0.681753551741409 0.824516206571175;0.0362543176279101 0.689080615156419 0.821454859237968;0.0323879761904762 0.696407678571429 0.818393511904762;0.028790380952381 0.700346142857143 0.81472280952381;0.0251927857142857 0.704284607142857 0.811052107142857;0.0215951904761905 0.708223071428571 0.807381404761905;0.0179975952380952 0.712161535714286 0.803710702380952;0.0144 0.7161 0.80004;0.0099 0.719 0.79386;0.0054 0.7219 0.78768;0.000900000000000001 0.7248 0.7815;0.0018 0.7275 0.7752;0.0046 0.7301 0.7688;0.0094 0.7327 0.7623;0.0162 0.7352 0.7558;0.0253 0.7376 0.7492;0.0369 0.74 0.7426;0.0504 0.7423 0.7359;0.0638 0.7446 0.7292;0.077 0.7468 0.7224;0.0899 0.7489 0.7156;0.1023 0.751 0.7088;0.1141 0.7531 0.7019;0.1252 0.7552 0.695;0.1354 0.7572 0.6881;0.1448 0.7593 0.6812;0.1532 0.7614 0.6741;0.1609 0.7635 0.6671;0.1678 0.7656 0.6599;0.1741 0.7678 0.6527;0.1799 0.7699 0.6454;0.1853 0.7721 0.6379;0.1905 0.7743 0.6303;0.1954 0.7765 0.6225;0.2003 0.7787 0.6146;0.2061 0.7808 0.6065;0.2118 0.7828 0.5983;0.2178 0.7849 0.5899;0.2244 0.7869 0.5813;0.2318 0.7887 0.5725;0.2401 0.7905 0.5636;0.2491 0.7922 0.5546;0.2589 0.7937 0.5454;0.2695 0.7951 0.536;0.2809 0.7964 0.5266;0.2929 0.7975 0.517;0.3052 0.7985 0.5074;0.3176 0.7994 0.4975;0.3301 0.8002 0.4876;0.3424 0.8009 0.4774;0.3548 0.8016 0.4669;0.3671 0.8021 0.4563;0.3795 0.8026 0.4454;0.3921 0.8029 0.4344;0.405 0.8031 0.4233;0.4184 0.803 0.4122;0.4322 0.8028 0.4013;0.4463 0.8024 0.3904;0.4608 0.8018 0.3797;0.4753 0.8011 0.3691;0.4899 0.8002 0.3586;0.5044 0.7993 0.348;0.5187 0.7982 0.3374;0.5329 0.797 0.3267;0.547 0.7957 0.3159;0.5609 0.7943 0.305;0.5748 0.7929 0.2941;0.5886 0.7913 0.2833;0.6024 0.7896 0.2726;0.6161 0.7878 0.2622;0.6297 0.7859 0.2521;0.6433 0.7839 0.2423;0.6567 0.7818 0.2329;0.6701 0.7796 0.2239;0.6833 0.7773 0.2155;0.6963 0.775 0.2075;0.7091 0.7727 0.1998;0.7218 0.7703 0.1924;0.7344 0.7679 0.1852;0.7468 0.7654 0.1782;0.759 0.7629 0.1717;0.771 0.7604 0.1658;0.7829 0.7579 0.1608;0.7945 0.7554 0.157;0.806 0.7529 0.1546;0.8172 0.7505 0.1535;0.8281 0.7481 0.1536;0.8389 0.7457 0.1546;0.8495 0.7435 0.1564;0.86 0.7413 0.1587;0.8703 0.7392 0.1615;0.8804 0.7372 0.165;0.8903 0.7353 0.1695;0.9 0.7336 0.1749;0.9093 0.7321 0.1815;0.9184 0.7308 0.189;0.9272 0.7298 0.1973;0.9357 0.729 0.2061;0.944 0.7285 0.2151;0.9523 0.7284 0.2237;0.9606 0.7285 0.2312;0.9689 0.7292 0.2373;0.977 0.7304 0.2418;0.9842 0.733 0.2446;0.99 0.7365 0.2429;0.9946 0.7407 0.2394;0.9966 0.7458 0.2351;0.9971 0.7513 0.2309;0.9972 0.7569 0.2267;0.9971 0.7626 0.2224;0.9969 0.7683 0.2181;0.9966 0.774 0.2138;0.9962 0.7798 0.2095;0.9957 0.7856 0.2053;0.9949 0.7915 0.2012;0.9938 0.7974 0.1974;0.9923 0.8034 0.1939;0.9906 0.8095 0.1906;0.9885 0.8156 0.1875;0.9861 0.8218 0.1846;0.9835 0.828 0.1817;0.9807 0.8342 0.1787;0.9778 0.8404 0.1757;0.9748 0.8467 0.1726;0.972 0.8529 0.1695;0.9694 0.8591 0.1665;0.9671 0.8654 0.1636;0.9651 0.8716 0.1608;0.9634 0.8778 0.1582;0.9619 0.884 0.1557;0.9608 0.8902 0.1532;0.9601 0.8963 0.1507;0.9596 0.9023 0.148;0.9595 0.9084 0.145;0.9597 0.9143 0.1418;0.9601 0.9203 0.1382;0.9608 0.9262 0.1344;0.9618 0.932 0.1304;0.9629 0.9379 0.1261;0.9642 0.9437 0.1216;0.9657 0.9494 0.1168;0.9674 0.9552 0.1116;0.9692 0.9609 0.1061;0.9711 0.9667 0.1001;0.973 0.9724 0.0938;0.9749 0.9782 0.0872;0.9769 0.9839 0.0805],...
%     'FontSize',24,'Layer','top','LineWidth',1.5);

% fprintf("binding energy of excitations at k=0 (in the unit of t):");
% disp(omega_set(islocalmax(A(:,1))) -mu)