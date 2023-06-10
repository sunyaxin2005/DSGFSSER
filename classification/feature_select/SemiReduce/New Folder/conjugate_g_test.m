
function [perish,p]=conjugate_g_test(train_mat,test_mat,R,p)
[~,n]=size(train_mat);
pcom=zeros(n,1);
xicom=zeros(n,1);
ftol = 0.0000001;
[p]=frprmn(R,train_mat,test_mat,pcom,xicom,p, ftol);
[~,perish]=functionTwo(train_mat,test_mat,R,p,1);
end
function a=sgn(x)
if x>0
    a=1;
end
if x==0;
    a=0;
end
if x<0
    a=-1;
end
end
function [sum,perish]=functionTwo(train_mat,test_mat,R,xt,a)
	sum1=0;
	sum2=0;

    m=size(xt);
    m1=size(test_mat);
    
	for i=1:1:m
		sum1=sum1+abs(xt(i));
    end
	
	%train_mat*weight_mat
	product=train_mat*xt;
	
	for i=1:1:m1
		I1=test_mat(i);
        I2=product(i);
		sum2=sum2+((I1-I2)^2);
    end
	
	if(a==1)
		perish=sum2;
    else
        if(a==0)
		    perish=1000;
        else
            disp('error in a');
        end
    end

	sum=sum1+R*sum2;
end
function [num]=functionOne(train_mat,test_mat,pcom,xicom,x,R)
	num=f1dim(train_mat,test_mat,pcom,xicom,x,R);
end
function [num]=f1dim(train_mat,test_mat,pcom,xicom,x,R)
	[m,n]=size(train_mat);
    [m1,n1]=size(test_mat);
    xt=zeros(n,1);
	for j=1:1:n
        xt(j)=pcom(j)+x*xicom(j);
    end
    
    [num,perish]=functionTwo(train_mat,test_mat,R,xt,0);
end

function [obstacle]=obstaclePara(train_mat,test_mat,p)	
	product=train_mat*p;
    [m,n]=size(test_mat);
    obstacle=zeros(m,1);
	for i=1:1:m
		I1=test_mat(i);
        I2=product(i);
        I_MINUS=I1-I2;%%%%%%%%%%%%%%%%%%%%%%%
		obstacle(i)=I_MINUS;
    end
end
function [gra]=gradient(train_mat,test_mat,p,R)
	obstacle=obstaclePara(train_mat,test_mat,p);  
    [m,n]=size(train_mat);
    gra=zeros(n,1);
    
	for i=1:1:n
		sum=0;
		for j=1:1:m
			I=train_mat(j,i);
		    k=obstacle(j);
			sum=sum+I*k;
        end
		if p(i)>=0
            gra(i)=1-2*R*sum;
        end
        if p(i)<0
            gra(i)=(-1)-2*R*sum;
        end
    end
end
function [fret,pcom,xicom,xi,p]=functionMin(train_mat,test_mat,xicom,p,pcom,xi,R)  
    tol = 0.000000001;
   [ncom,mcom] =size(p);
    for j=1:1:ncom
		pcom(j)=p(j);
        xicom(j)=xi(j);
    end
    ax = 0.0;
    xx = 1.0;
    bx=1;
    fa=1;
    fb=1;
    fx=1;
    xmin=1;
    [ax,xx,bx,fa,fx,fb]=mnbrak(train_mat,test_mat,pcom,xicom,ax,xx,bx,fa,fx,fb,R);
    [xmin,fret] = brent(train_mat,test_mat,pcom,xicom,ax, xx, bx, tol,xmin,R);
    for j =1:1:ncom
        xi(j)=xmin*xi(j);
        p(j)=p(j)+xi(j);
    end
end
function [ax,bx,cx,fa,fb,fc]=mnbrak(train_mat,test_mat,pcom,xicom,ax,bx,cx,fa,fb,fc,R)
    gold = 1.618034;
    glimit = 100;
	tiny = 1e-20;
    fa = functionOne(train_mat,test_mat,pcom,xicom,ax,R);
    fb = functionOne(train_mat,test_mat,pcom,xicom,bx,R);
    if fb > fa
        dum = ax;
        ax = bx;
        bx = dum;
        dum = fb;
        fb = fa;
        fa = dum;
    end
    cx = bx + gold * (bx - ax);
    fc = functionOne(train_mat,test_mat,pcom,xicom,cx,R);
	while fb >= fc
        r = (bx - ax) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        dum = q - r;
        if abs(dum) < tiny
			dum = tiny;
        end
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2 * dum);
        ulim = bx + glimit * (cx - bx);
        if (bx - u) * (u - cx) > 0
            fu = functionOne(train_mat,test_mat,pcom,xicom,u,R);
            if fu < fc
                ax = bx;
                fa = fb;
                bx = u;
                fb = fu;
                break;
            else
				if fu > fb
					cx = u;
					fc = fu;
					break;
                end
            end    
    
            u = cx + gold * (cx - bx);
            fu = functionOne(train_mat,test_mat,pcom,xicom,u,R);
        else
			if (cx - u) * (u - ulim) > 0
				fu = functionOne(train_mat,test_mat,pcom,xicom,u,R);
				if fu < fc
					bx = cx;
					cx = u;
					u = cx + gold * (cx - bx);
					fb = fc;
					fc = fu;
					fu = functionOne(train_mat,test_mat,pcom,xicom,u,R);
                end
            else
				if (u - ulim) * (ulim - cx) >= 0
					u = ulim;
					fu = functionOne(train_mat,test_mat,pcom,xicom,u,R);
                else
					u = cx + gold * (cx - bx);
					fu = functionOne(train_mat,test_mat,pcom,xicom,u,R);
                end
            end
        end
		ax = bx;
		bx = cx;
		cx = u;
		fa = fb;
		fb = fc;
		fc = fu;
    end
end
function [xmin,fret]=brent(train_mat,test_mat,pcom,xicom,ax,bx,cx,tol,xmin,R)
    itmax = 100;
	cgold = 0.381966;
	zeps = 0.000000001;
    a = ax;
    if cx < ax
		a = cx;
    end
    b = ax;
    if cx > ax
		b = cx;
    end
    v = bx;
    w = v;
    x = v;
    e = 0.0;
    fx = functionOne(train_mat,test_mat,pcom,xicom,x,R);
    fv1 = fx;
    fw = fx;
    for iter =2:1:(itmax+1)
        xm = 0.5 * (a + b);
        tol1 = tol * abs(x) + zeps;
        tol2 = 2.0 * tol1;
        if abs(x - xm) <= tol2 - 0.5 * (b - a)
			break;
        end
        done = -1;
        if abs(e) > tol1
            r = (x - w) * (fx - fv1);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0
				p = -p;
            end
            q = abs(q);
            etemp = e;
            e = d;
            dum = abs(0.5 * q * etemp);
            if abs(p) < dum && p > q * (a - x) && p < q * (b - x)
                d = p / q;
                u = x + d;
                if u - a < tol2 || b - u < tol2
                    d = abs(tol1) * sgn(xm - x);
                end
                done = 0;
            end
        end
        if done
            if x >= xm
                e = a - x;
            else
                e = b - x;
            end
            d = cgold * e;
        end
        if abs(d) >= tol1
            u = x + d;
        else
            u = x + abs(tol1) * sgn(d);
        end
        fu = functionOne(train_mat,test_mat,pcom,xicom,u,R);
        if fu <= fx
            if u >= x
                a = x;
            else
                b = x;
            end
            v = w;
            fv1 = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        else
            if u < x
                a = u;
            else
                b = u;
            end
            if fu <= fw || w == x
                v = w;
                fv1 = fw;
                w = u;
                fw = fu;
            else
				if fu <= fv1 || v == x || v == w
					v = u;
					fv1 = fu;
                end
            end
        end
        
    end
    xmin = x;
    fret=fx;
end
function [p]=frprmn(R,train_mat,test_mat,pcom,xicom,p,ftol)
	itmax = 100;
    eps = 0.000000000001;
    [~,n1]=size(train_mat);
    mm=size(p);
    g=zeros(n1,1);
    h=zeros(n1,1);
    xi=zeros(n1,1);
	
	[fp,perish] = functionTwo(train_mat,test_mat,R,p,0);
    [xi]=gradient(train_mat,test_mat,p,R);
	for j=1:1:mm
        g(j)=(-1)*xi(j);
        h(j)=g(j);
        xi(j)=h(j);
    end
    for its=1:1:itmax
		if abs(perish)<0.00001
			break;
        end
        [fret,pcom,xicom,xi,p]=functionMin(train_mat,test_mat,xicom,p,pcom,xi,R);
		
        if 2.0 * abs(fret - fp) <= ftol * (abs(fret) + abs(fp) +eps)
            break;
        end 
        [fp,perish] = functionTwo(train_mat,test_mat,R,p,0);
        [xi]=gradient(train_mat,test_mat,p,R);
        gg = 0.0;
        dgg = 0.0;
		for j=1:1:mm
            gg = gg + g(j)*g(j);
            %dgg = dgg + CV_MAT_ELEM(*xi,double,j,0)*CV_MAT_ELEM(*xi,double,j,0);        // polak-ribiere ·¨
            dgg = dgg + (xi(j)+  g(j)) *  xi(j); 
        end
        if gg == 0.0
			break;
        end
         gam = dgg / gg;
        for j=1:1:mm
            g(j)=(-1)*xi(j);
            h(j)=g(j)+gam*h(j);
            xi(j)=h(j);
        end
    end	
end




