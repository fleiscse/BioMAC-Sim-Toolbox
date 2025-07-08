function speedR  = plotBeltSpeed(obj, X)

    bR = obj.idx.belt_right;
    speedR = X(bR);
    plot(speedR)

end