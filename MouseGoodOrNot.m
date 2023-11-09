function MouseGoodOrNot(src, e)



switch e.Source.SelectionType
    case 'normal'
        src.UserData = 1;
        disp('Good data');
    case 'alt'
        src.UserData = 2;
        disp('Bad data');
end
end

