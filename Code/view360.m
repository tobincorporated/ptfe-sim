t=1;


for th=180:-2:180-360;
    pause(1);
    view(th,0);
    time=num2str(t);
    saveas(gcf, ['animation/animation' time],'jpg');
    t=t+1;
end
% 
% for el=20:380
%     pause(.05);
%     view(20,el);
%     time=num2str(t);
%     saveas(gcf, ['animation/animation' time],'jpg');
%     t=t+1;
% end