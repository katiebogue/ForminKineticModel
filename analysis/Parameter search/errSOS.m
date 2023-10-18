function SOS=errSOS(model,exp,errtop,errbot)
    SOS=0;
    for i=1:length(model)
        max=exp(i)+errtop(i);
        min=exp(i)-errbot(i);
        if (model(i)<=max) && (model(i)>=min)
            dif=0;
        else
            dif=abs((model(i)-exp(i))^2);
        end
        SOS= SOS + dif;
    end
end