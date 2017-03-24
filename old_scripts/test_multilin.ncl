       y  = (/40,50,50,70,65,65,80/)
       N  = dimsizes(y)
     
       M  = 2
       x1 = (/100,200,300,400,500,600,700/)
       x2 = (/ 10, 20, 10, 30, 20, 20, 30/)
                                              ; create independent array
       X  = new ( (/M+1,N/), "float" )  
       X(0,:) = 1.0                           ; constant term on rhs
       X(1,:) = x1
       X(2,:) = x2
                                              ; partial regression coef
       beta = reg_multlin (y,X,False)
       print(beta)
