function [imgs] = generate_large_imgs(hea_vac,Pars)

newimgs = zeros(Pars.pixel_Nx,Pars.pixel_Nx,1,length(hea_vac(:,1)));

for i = 1:length(hea_vac(:,1))
    
    newimgs(11:14,11:14,1,i)= ones(4)*hea_vac(i,55); % Vacancy/Defect
    
    newimgs(9:10,  9:10,1,i)= ones(2)*hea_vac(i,1);
    newimgs(9:10, 11:12,1,i)= ones(2)*hea_vac(i,2);
    newimgs(9:10, 13:14,1,i)= ones(2)*hea_vac(i,3);
    newimgs(9:10, 15:16,1,i)= ones(2)*hea_vac(i,4);

    newimgs(11:12, 15:16,1,i)= ones(2)*hea_vac(i,5);
    newimgs(13:14, 15:16,1,i)= ones(2)*hea_vac(i,6);
    newimgs(15:16, 15:16,1,i)= ones(2)*hea_vac(i,7);

    newimgs(15:16, 13:14,1,i)= ones(2)*hea_vac(i,8);
    newimgs(15:16, 11:12,1,i)= ones(2)*hea_vac(i,9);
    newimgs(15:16,  9:10,1,i)= ones(2)*hea_vac(i,10);

    newimgs(13:14,  9:10,1,i)= ones(2)*hea_vac(i,11);
    newimgs( 9:10, 11:12,1,i)= ones(2)*hea_vac(i,12);

    newimgs( 7:8,  11:14,1,i)= ones(2,4)*hea_vac(i,13);
    newimgs(11:14, 17:18,1,i)= ones(4,2)*hea_vac(i,14);
    newimgs(17:18, 11:14,1,i)= ones(2,4)*hea_vac(i,15);
    newimgs(11:14,  7:8 ,1,i)= ones(4,2)*hea_vac(i,16);

    newimgs( 7:8,   7:8,1,i)= ones(2,2)*hea_vac(i,17);
    newimgs(17:18,17:18,1,i)= ones(2,2)*hea_vac(i,17);

    newimgs( 7:8, 17:18,1,i)= ones(2,2)*hea_vac(i,18);
    newimgs(17:18, 7:8 ,1,i)= ones(2,2)*hea_vac(i,18);

    newimgs( 3:6, 7:8 ,1,i)= ones(4,2)*hea_vac(i,19);
    newimgs( 3:6, 9:10,1,i)= ones(4,2)*hea_vac(i,20);
    newimgs( 3:6,11:12,1,i)= ones(4,2)*hea_vac(i,21);
    newimgs( 3:6,13:14,1,i)= ones(4,2)*hea_vac(i,22);
    newimgs( 3:6,15:16,1,i)= ones(4,2)*hea_vac(i,23);
    newimgs( 3:6,17:18,1,i)= ones(4,2)*hea_vac(i,24);

    newimgs( 7:8,19:22,1,i)= ones(2,4)*hea_vac(i,25);
    newimgs( 9:10,19:22,1,i)= ones(2,4)*hea_vac(i,26);
    newimgs(11:12,19:22,1,i)= ones(2,4)*hea_vac(i,27);
    newimgs(13:14,19:22,1,i)= ones(2,4)*hea_vac(i,28);
    newimgs(15:16,19:22,1,i)= ones(2,4)*hea_vac(i,29);
    newimgs(17:18,19:22,1,i)= ones(2,4)*hea_vac(i,30);

    newimgs(19:22,17:18 ,1,i)= ones(4,2)*hea_vac(i,31);
    newimgs(19:22,15:16 ,1,i)= ones(4,2)*hea_vac(i,32);
    newimgs(19:22,13:14 ,1,i)= ones(4,2)*hea_vac(i,33);
    newimgs(19:22,11:12 ,1,i)= ones(4,2)*hea_vac(i,34);
    newimgs(19:22, 9:10 ,1,i)= ones(4,2)*hea_vac(i,35);
    newimgs(19:22, 7:8  ,1,i)= ones(4,2)*hea_vac(i,36);

    newimgs(17:18,3:6,1,i)= ones(2,4)*hea_vac(i,37);
    newimgs(15:16,3:6,1,i)= ones(2,4)*hea_vac(i,38);
    newimgs(13:14,3:6,1,i)= ones(2,4)*hea_vac(i,39);
    newimgs(11:12,3:6,1,i)= ones(2,4)*hea_vac(i,40);
    newimgs( 9:10,3:6,1,i)= ones(2,4)*hea_vac(i,41);
    newimgs( 7:8,3:6,1,i) = ones(2,4)*hea_vac(i,42);

    newimgs( 1:2,7:12,1,i) = ones(2,6)*hea_vac(i,43);
    newimgs( 1:2,13:18,1,i) = ones(2,6)*hea_vac(i,44);

    newimgs( 7:12,23:24,1,i) = ones(6,2)*hea_vac(i,45);
    newimgs(13:18,23:24,1,i) = ones(6,2)*hea_vac(i,46);

    newimgs(23:24,13:18,1,i) = ones(2,6)*hea_vac(i,47);
    newimgs(23:24, 7:12,1,i) = ones(2,6)*hea_vac(i,48);

    newimgs(13:18,1:2,1,i) = ones(6,2)*hea_vac(i,49);
    newimgs( 7:12,1:2,1,i) = ones(6,2)*hea_vac(i,50);

    newimgs( 1:4,1:2,1,i) = ones(4,2)*hea_vac(i,51);
    newimgs( 1:2,3:4,1,i) = ones(2,2)*hea_vac(i,51);

    newimgs( 1:4,23:24,1,i) = ones(4,2)*hea_vac(i,52);
    newimgs( 1:2,21:22,1,i) = ones(2,2)*hea_vac(i,52);

    newimgs(21:24,23:24,1,i) = ones(4,2)*hea_vac(i,53);
    newimgs(23:24,21:22,1,i) = ones(2,2)*hea_vac(i,53);

    newimgs(21:24,1:2,1,i) = ones(4,2)*hea_vac(i,54);
    newimgs(23:24,3:4,1,i) = ones(2,2)*hea_vac(i,54);
    
    
end

imgs = newimgs;

end

