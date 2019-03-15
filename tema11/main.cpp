#include <iostream>
using namespace std;

class NumarComplex
{
private:
    double re,im;
public:
    NumarComplex(double real=0,double imag=0)
    {
        this->im=imag;
        this->re=real;
    }
    friend istream& operator >>(istream& f,NumarComplex& z);
    friend ostream& operator <<(ostream& g,NumarComplex& z);
    NumarComplex& operator+(NumarComplex& z);
    NumarComplex& operator*(NumarComplex& z);
    NumarComplex& operator/(NumarComplex& z);
    void setReal(double real)
    {
        this->re=real;
    }
    double getReal()
    {
        return this->re;
    }
    double getImag()
    {
        return this->im;
    }
    friend class Matrice;
};
istream& operator >>(istream& f,NumarComplex& z)
{
    cout<<"Introduceti partea reala: ";
    f>>z.re;
    cout<<"Introduceti partea imaginara: ";
    f>>z.im;
    return f;
}
ostream& operator<<(ostream& g,NumarComplex& z)
{
    g<<"z= "<<z.re<<"+"<<z.im<<"i"<<" ";
    return g;
}
NumarComplex& NumarComplex::operator+(NumarComplex& z)
{
    NumarComplex *z1=new NumarComplex();
    z1->re= this->re+ z.re;
    z1->im= this->im+ z.im;
    return *z1;
}
NumarComplex& NumarComplex::operator*(NumarComplex& z)
{
    NumarComplex *z1=new NumarComplex();
    z1->re= (this->re* z.re-this->im*z.im);
    z1->im= (this->re*z.im+this->im*z.re);
    return *z1;
}
NumarComplex& NumarComplex::operator/(NumarComplex& z)
{
    NumarComplex *z1=new NumarComplex();
    z1->re=(this->re*z.re+this->im*z.im)/(z.re*z.re+z.im*z.im);
    z1->im=(z.re*this->im-this->re*z.im)/(z.re*z.re+z.im*z.im);
    return *z1;
}
class Matrice
{
private:
    int n,m;
    NumarComplex **v;
    void detcofactor(NumarComplex **mat,NumarComplex **&temp,int p,int q,int n);
    void adjuncta(NumarComplex **&mat);
public:
    Matrice(int lin,int col)
    {
        this->n=lin;
        this->m=col;
        v=new NumarComplex*[n]();
        for(int i=0; i<n; i++)
            v[i]=new NumarComplex[m]();
    }
    NumarComplex** getv()
    {
        return this->v;
    }
    int getn()
    {
        return this->n;
    }
    friend istream& operator>>(istream& f,Matrice& mat);
    friend ostream& operator<<(ostream& g,Matrice& mat);
    Matrice& operator+(Matrice& mat);
    Matrice& operator*(Matrice& mat);
    NumarComplex& determinant(NumarComplex **mat,int n);
    Matrice& inversa();
};
istream& operator>>(istream& f,Matrice& mat)
{
    cout<< "Introduceti matricea:";
    for(int i=0; i<mat.n; i++)
    {
        for(int j=0; j<mat.m; j++)
            f>>mat.v[i][j];
    }
    return f;
}
ostream& operator<<(ostream& g,Matrice& mat)
{
    for(int i=0; i<mat.n; i++)
    {
        for(int j=0; j<mat.m; j++)
            g<<mat.v[i][j];
        g<<endl;
    }
    return g;
}
Matrice& Matrice::operator+(Matrice& mat)
{
    Matrice *m=new Matrice(mat.n,mat.m);
    for(int i=0; i<mat.n; i++)
    {
        for(int j=0; j<mat.m; j++)
            m->v[i][j]=this->v[i][j]+mat.v[i][j];
    }
    return *m;
}
Matrice& Matrice::operator*(Matrice& mat)
{
    Matrice *m=new Matrice(mat.n,mat.m);
    for(int i=0; i<this->m; i++)
    {
        for(int j=0; j<this->m; j++)
            for(int k=0; k<this->m; k++)
                m->v[i][j]=m->v[i][j]+this->v[i][k]*mat.v[k][j];
    }
    return *m;
}
void Matrice::detcofactor(NumarComplex **mat,NumarComplex **&temp,int p,int q,int n)
{
    int i=0;
    int j=0;
    temp=new NumarComplex *[n];
    for(int k=0; k<n; k++)
        temp[k]=new NumarComplex[n]();
    for(int linie=0; linie<n; linie++)
    {
        for(int col=0; col<n; col++)
        {
            if(linie!=p && col!=q)
            {
                temp[i][j++]=mat[linie][col];
                if(j==n-1)
                {
                    j=0;
                    i++;
                }
            }
        }
    }
}
NumarComplex& Matrice::determinant(NumarComplex **mat,int n)
{
    NumarComplex *rezultat= new NumarComplex;
    if(n==1)
        return mat[0][0];
    NumarComplex **temp;
    int sgn=1;
    NumarComplex sign;
    sign.setReal(sgn);
    for(int i=0; i<n; i++)
    {
        this->detcofactor(mat,temp,0,i,n);
        *rezultat=*rezultat+sign*mat[0][i]*determinant(temp,n-1);
        sgn=-sgn;
        sign.setReal(sgn);
    }
    return *rezultat;

}
void Matrice::adjuncta(NumarComplex **&mat)
{
    if(this->getn()==1)
    {
        mat[0][0]=1;
        return;
    }
    NumarComplex sign(1,0), **temp;
    int n=this->getn();
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            this->detcofactor(this->getv(),temp,i,j,n);
            if((i+j)%2==0)
                sign.setReal(1);
            else
                sign.setReal(-1);
            mat[j][i]=sign*this->determinant(temp,n-1);
        }
    }
}
Matrice& Matrice::inversa()
{
    NumarComplex det=this->determinant(this->getv(),this->getn());
    if(det.getReal()==0 && det.getImag()==0)
    {
        cout<<"Determinantul este 0"<<endl;
        return *this;
    }
    int n=this->getn();
    Matrice *m_inversa=new Matrice(n,n);
    NumarComplex **adj;
    adj=new NumarComplex*[n];
    for(int i=0; i<this->getn(); i++)
        adj[i]=new NumarComplex[n]();
    this->adjuncta(adj);
    for(int i=0; i<n; i++)
        for(int j=0; j<n; j++)
            m_inversa->v[i][j]=adj[i][j]/det;
return *m_inversa;
        }
int main()
{
//NumarComplex z1(2,3),z2;
//cin>> z2;
//cout<<z1<<endl;
//cout<<z2<<endl;
//cout<<z1+z2<<endl;
//cout<<z1*z2<<endl;

Matrice m1(2,2),m2(2,2),m3(2,3);
cin>>m1; cin>>m3;
cout<<m1<<endl<<m3<<endl;
//cout<<m1+m2<<endl;
cout<<m1*m3<<endl;
//cout<<m1.determinant(m1.getv(),m1.getn())<<endl;
//cout<<m1.inversa();

    return 0;
}
