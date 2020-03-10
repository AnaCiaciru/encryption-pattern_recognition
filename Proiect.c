#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct pixel {
    unsigned char r, g, b;
};

struct detectii {
    unsigned int fi, fj;
    double corre;
    struct pixel C;
};


void incarca_imagine(char *caleImagine, unsigned int *H, unsigned int *W, unsigned char ***v, unsigned char **header) {
    FILE *in = fopen(caleImagine, "rb");
    if (in == NULL) {
        printf("Nu s-a deschis fisierul de intrare");
        exit(0);
    }

    ///pastram header-ul
    *header = (unsigned char *) malloc(54 * sizeof(unsigned char));

    unsigned int k, dimOct;

    for (k = 0; k < 54; k++)
        fread(*header + k, sizeof(unsigned char), 1, in);

    fseek(in, 2, SEEK_SET);
    fread(&dimOct, sizeof(unsigned int), 1, in);

    /// W x H in pixeli
    fseek(in, 18, SEEK_SET);
    fread(W, sizeof(unsigned int), 1, in);

    fseek(in, 22, SEEK_SET);
    fread(H, sizeof(unsigned int), 1, in);

    unsigned int padd = 0, dw, dv;
    int i, j;

    dw = 3 * (*W);///dimensiunea in octeti a unei linii
    if (dw % 4 > 0) padd = 4 - dw % 4;

    fseek(in, 54, SEEK_SET); ///trecem peste header

    ///alocam memorie
    dv = (*H) * (*W); ///dimensiunea matricei liniarizate
    *v = (unsigned char **) malloc(dv * sizeof(unsigned char *));
    for (i = 0; i < dv; i++)
        *(*v + i) = (unsigned char *) malloc(3 * sizeof(unsigned char)); ///pentru fiecare pixel
    unsigned char x;///pentru a citi octetii de padding
    ///citim in vector
    for (i = *H - 1; i >= 0; i--) {
        for (j = 0; j < *W; j++)
            for (k = 0; k < 3; k++)
                fread(*(*v + (*W) * i + j) + k, sizeof(unsigned char), 1, in);
        for (k = 0; k < padd; k++)
            fread(&x, sizeof(unsigned char), 1, in);
    }
    fclose(in);
}

void descarca_imagine(char *caleImgExtern, unsigned int H, unsigned int W, unsigned char **v, unsigned char *header) {
    FILE *out = fopen(caleImgExtern, "wb");
    if (out == NULL) {
        printf("Nu s-a deschis fisierul de iesire");
        exit(0);
    }

    unsigned int k, padd = 0, dw;
    int i, j;
    for (k = 0; k < 54; k++) {
        fwrite(header + k, sizeof(unsigned char), 1, out);
        fflush(out);
    }

    dw = 3 * W; ///dimensiunea in octeti a unei linii
    if (dw % 4 > 0) padd = 4 - dw % 4;
    unsigned char zero;
    zero = '0';
    for (i = H - 1; i >= 0; i--) {
        for (j = 0; j < W; j++)
            for (k = 0; k < 3; k++) {
                fwrite(*(v + W * i + j) + k, sizeof(unsigned char), 1, out);
                fflush(out);
            }
        for (k = 0; k < padd; k++) {
            fwrite(&zero, sizeof(unsigned char), 1, out);
            fflush(out);
        }
    }
    fclose(out);
}

unsigned int xorshift(unsigned int *seed) {
    unsigned int x = *seed;
    x = x ^ (x << 13);
    x = x ^ (x >> 17);
    x = x ^ (x << 5);
    *seed = x;
    return x;
}

void xorshift32(unsigned int seed, unsigned int **v, unsigned int n) {
    *v = (unsigned int *) malloc(n * sizeof(unsigned int));
    unsigned int i;
    for (i = 1; i <= n; i++)
        *(*v + i) = xorshift(&seed);
}

void permutare(unsigned int n, unsigned int **perm, unsigned int R0, unsigned int **rand) {
    unsigned int m, r, aux, i, k;
    ///n= W*H
    m = 2 * n - 1;
    xorshift32(R0, rand, m);

    *perm = (unsigned int *) malloc(n * sizeof(unsigned int));

    for (i = 0; i < n; i++) *(*perm + i) = i;
    k = 1;
    for (i = n - 1; i >= 1; i--) {
        r = (*(*rand + k)) % (i + 1);
        k++;
        aux = *(*perm + i);
        *(*perm + i) = *(*perm + r);
        *(*perm + r) = aux;
    }
}

void permutare_imagine(unsigned int n, unsigned int *perm, unsigned char ***v) {
    unsigned int i, j, k;
    unsigned char **vv;
    ///aloc memorie
    vv = (unsigned char **) malloc(n * sizeof(unsigned char *));
    for (i = 0; i < n; i++)
        vv[i] = (unsigned char *) malloc(3 * sizeof(unsigned char)); ///pentru fiecare pixel

    ///modificam pixelii conform permutarii
    for (i = 0; i < n; i++) {
        k = perm[i];
        for (j = 0; j < 3; j++)
            vv[k][j] = (*v)[i][j];
    }
    ///copiem in vectorul initial rezultatul
    for (i = 0; i < n; i++)
        for (j = 0; j < 3; j++)
            (*v)[i][j] = vv[i][j];
    ///eliberez memoria utilizata auxiliar in vv
    for (i = 0; i < n; i++) free(*(vv + i));
    free(vv);
}

unsigned char* xor(
unsigned char *p1,
unsigned char *p2
)
{
unsigned int i;
unsigned char *p;
p = (unsigned char *) malloc(3 * sizeof(unsigned char));
for (
i = 0;
i<3;i++)
p[i]= p1[i] ^ p2[i];
return
p;
}

unsigned char *xor_int(unsigned char *p, unsigned int x) {
    unsigned char *t, *q;
    int i;
    q = (unsigned char *) malloc(3 * sizeof(unsigned char));
    t = (&x);
    for (i = 0; i < 3; i++) {
        q[i] = p[i] ^ (*t);
        t = t + 1;
    }
    return q;
}

void xorare(unsigned char ***v, unsigned int n, unsigned int sv, unsigned int *rand) {
    unsigned int i;
    (*v)[0] = xor_int(xor_int((*v)[0], sv), rand[n]);///cazul imediat
    for (i = 1; i < n; i++)
        (*v)[i] = xor_int(xor((*v)[i - 1], (*v)[i]), rand[n + i] ); ///recurenta
}


void criptare(char *caleImagine, char *caleImgExtern, char *keyText) {
    unsigned int H, W, n, *perm, i, j, key, *rand, sv;
    unsigned char **v, *header;
    incarca_imagine(caleImagine, &H, &W, &v, &header);
    n = H * W;

    FILE *in = fopen(keyText, "r");
    if (in == NULL) {
        printf("Eroare la deschiderea cheii secrete");
        exit(0);
    }
    fscanf(in, "%u%u", &key, &sv);
    fclose(in);

    permutare(n, &perm, key, &rand);
    permutare_imagine(n, perm, &v);
    xorare(&v, n, sv, rand);
    descarca_imagine(caleImgExtern, H, W, v, header);
}

void permutare_inversa(unsigned int n, unsigned int **perm) {
    unsigned int *p, i;
    p = (unsigned int *) malloc(n * sizeof(unsigned int));
    for (i = 0; i < n; i++)
        p[(*perm)[i]] = i;
    for (i = 0; i < n; i++)
        (*perm)[i] = p[i];
    free(p);
}

void xorare1(unsigned char ***v, unsigned int n, unsigned int sv, unsigned int *rand) {
    unsigned int i, j;
    unsigned char **a;
    a = (unsigned char **) malloc(n * sizeof(unsigned char *));
    for (i = 0; i < n; i++)
        a[i] = (unsigned char *) malloc(3 * sizeof(unsigned char)); ///pentru fieca
    a[0] = xor_int(xor_int((*v)[0], sv), rand[n]);///cazul imediat
    for (i = 1; i < n; i++)
        a[i] = xor_int(xor((*v)[i - 1], (*v)[i]), rand[n + i] ); ///recurenta
    for (i = 0; i < n; i++)
        for (j = 0; j < 3; j++)
            (*v)[i][j] = a[i][j];
    ///eliberez memoria utilizata auxiliar in vv
    for (i = 0; i < n; i++) free(*(a + i));
    free(a);

}

void decriptare(char *cale_img_initiala, char *cale_img_criptata, char *keyText) {
    unsigned int H, W, key, sv, n, *perm, *rand;
    unsigned char **v, *header;
    incarca_imagine(cale_img_criptata, &H, &W, &v, &header);

    n = H * W;
    FILE *in = fopen(keyText, "r");
    if (in == NULL) {
        printf("Eroare la deschiderea cheii secrete");
        exit(0);
    }
    fscanf(in, "%u%u", &key, &sv);
    fclose(in);

    permutare(n, &perm, key, &rand);
    permutare_inversa(n, &perm);
    xorare1(&v, n, sv, rand);
    permutare_imagine(n, perm, &v);
    descarca_imagine(cale_img_initiala, H, W, v, header);
}


double calcul(unsigned int *f, double q) {
    int i;
    double x, sum = 0;
    for (i = 0; i < 256; i++) {
        x = f[i] - q;
        x = x * x;
        x /= q;
        sum += x;
    }
    return sum;
}

void chi__patrat(char *cale_imagine) {
    FILE *in = fopen(cale_imagine, "rb");
    if (in == NULL) {
        printf("Nu s-a deschis fisierul de intrare");
        exit(0);
    }
    unsigned int H, W;
    unsigned int *R, *G, *B;

    R = (unsigned int *) calloc(256, sizeof(unsigned int));
    G = (unsigned int *) calloc(256, sizeof(unsigned int));
    B = (unsigned int *) calloc(256, sizeof(unsigned int));

    /// W x H in pixeli
    fseek(in, 18, SEEK_SET);
    fread(&W, sizeof(unsigned int), 1, in);

    fseek(in, 22, SEEK_SET);
    fread(&H, sizeof(unsigned int), 1, in);

    unsigned int padd = 0, dw, dv;
    int i, j, k;

    dw = 3 * W;///dimensiunea in octeti a unei linii
    if (dw % 4 > 0) padd = 4 - dw % 4;

    fseek(in, 54, SEEK_SET); ///trecem peste header
    unsigned char x;///pentru a citi octetii de padding

    for (i = 0; i < H; i++) {
        for (j = 0; j < W; j++) {
            fread(&x, sizeof(unsigned char), 1, in);
            B[x]++;
            fread(&x, sizeof(unsigned char), 1, in);
            G[x]++;
            fread(&x, sizeof(unsigned char), 1, in);
            R[x]++;
        }
        for (k = 0; k < padd; k++)
            fread(&x, sizeof(unsigned char), 1, in);
    }
    fclose(in);
    double r, g, b;
    double q;
    q = W * H;
    q = q * (1.00) / 256;
    r = calcul(R, q);
    g = calcul(G, q);
    b = calcul(B, q);
    printf("%.2f\n%.2f\n%.2f\n", r, g, b);
}


void grayscale(char *cale_sursa, char *cale_destinatie) {
    FILE *in = fopen(cale_sursa, "rb");
    FILE *out = fopen(cale_destinatie, "wb");
    if (in == NULL) {
        printf("Fisierul sursa nu s-a putut deschide");
        exit(0);
    }
    unsigned int H, W, i, j;
    unsigned char R, G, B, zero;
    ///header
    unsigned char c;
    for (i = 0; i < 54; i++) {
        fread(&c, sizeof(unsigned char), 1, in);
        fwrite(&c, sizeof(unsigned char), 1, out);
    }
    fseek(in, 18, SEEK_SET);
    fread(&W, sizeof(unsigned int), 1, in);
    fread(&H, sizeof(unsigned int), 1, in);

    fseek(in, 54, SEEK_SET);

    ///padding

    unsigned int padd = 0;
    if ((3 * W) % 4 > 0)
        padd = 4 - (3 * W) % 4;

    ///body
    zero = '0';
    unsigned int x;
    for (i = 0; i < H; i++) {
        for (j = 0; j < W; j++) {
            fread(&B, sizeof(unsigned char), 1, in);
            fread(&G, sizeof(unsigned char), 1, in);
            fread(&R, sizeof(unsigned char), 1, in);
            x = (unsigned int) (0.299 * R + 0.587 * G + 0.114 * B);
            B = G = R = (unsigned char) x;
            fwrite(&B, sizeof(unsigned char), 3, out);
            fflush(out);
        }
        fseek(in, padd, SEEK_CUR);
        fwrite(&zero, sizeof(unsigned char), padd, out);
        fflush(out);
    }
    fclose(in);
    fclose(out);
}


void color(struct pixel **a, unsigned int h, unsigned int w, unsigned int fi, unsigned int fj, struct pixel C) {
    unsigned int i, j;
    for (i = fi; i < fi + h; i++)
        a[i][fj] = a[i][fj + w - 1] = C;

    for (j = fj; j < fj + w; j++)
        a[fi][j] = a[fi + h - 1][j] = C;
}

void colorare(char *imagine, unsigned int h, unsigned int w, unsigned int nv, struct detectii *v) {
    FILE *in = fopen(imagine, "rb+");
    if (in == NULL) {
        printf("Nu s-a deschis imaginea I");
        exit(0);
    }
    int t;
    unsigned int W, H, j, i;

    ///preluam datele din imaginea I
    fseek(in, 18, SEEK_SET);
    fread(&W, sizeof(unsigned int), 1, in);
    fread(&H, sizeof(unsigned int), 1, in);
    fseek(in, 54, SEEK_SET);

    ///padding I
    unsigned int padd = 0;
    if ((3 * W) % 4 > 0)
        padd = 4 - (3 * W) % 4;

    ///alocare matrice ptr I
    unsigned char zero = '0';
    struct pixel **a;
    a = (struct pixel **) malloc(H * sizeof(struct pixel *));
    for (i = 0; i < H; i++)
        *(a + i) = (struct pixel *) malloc(W * sizeof(struct pixel));

    ///citire matrice
    for (t = H - 1; t >= 0; t--) {
        for (j = 0; j < W; j++)
            fread(&(a[t][j]), sizeof(struct pixel), 1, in);
        fseek(in, padd, SEEK_CUR);
    }

    fseek(in, 54, SEEK_SET);

    for (i = 0; i < nv; i++)
        color(a, h, w, v[i].fi, v[i].fj, v[i].C);

    for (t = H - 1; t >= 0; t--) {
        for (j = 0; j < W; j++) {
            fwrite(&(a[t][j]), sizeof(struct pixel), 1, in);
            fflush(in);
        }
        fwrite(&zero, sizeof(unsigned char), padd, in);
        fflush(in);
    }

    for (i = 0; i < H; i++) free(a[i]);
    free(a);
    fclose(in);
}

void adaugare(unsigned int *n, struct detectii **v, unsigned int i, unsigned int j, double corr, struct pixel C) {
    struct detectii *aux;
    aux = (struct detectii *) realloc(*v, (*n + 1) * sizeof(struct detectii));
    if (aux == NULL) {
        printf("NU e suficienta memorie");
        free(*v);
        exit(0);
    } else {
        (*v) = aux;
        (*(*v + (*n))).fi = i;
        (*(*v + (*n))).fj = j;
        (*(*v + (*n))).corre = corr;
        (*(*v + (*n))).C = C;
        (*n)++;
    }
}


void template_matching(char *cale_imagine, char *cale_sablon, double ps, unsigned *nv, struct detectii **v,
                       struct pixel C) {   ///transformam imaginea si sablonul in tonuri de gri
    char grayImg[] = "gray_imagine.bmp";
    char graySablon[] = "gray_sablon.bmp";

    grayscale(cale_imagine, grayImg);
    grayscale(cale_sablon, graySablon);

    ///citim datele despre imagine
    FILE *in = fopen(grayImg, "rb");

    if (in == NULL) {
        printf("Nu s-a deschis imaginea I");
        exit(0);
    }
    int t;
    unsigned int W, H, j, w, h, i;
    ///preluam datele din imaginea I
    fseek(in, 18, SEEK_SET);
    fread(&W, sizeof(unsigned int), 1, in);
    fread(&H, sizeof(unsigned int), 1, in);
    fseek(in, 54, SEEK_SET);


    ///padding I
    unsigned int padd = 0;
    if ((3 * W) % 4 > 0)
        padd = 4 - (3 * W) % 4;

    ///alocare matrice ptr I
    unsigned char **a, x;
    a = (unsigned char **) malloc(H * sizeof(unsigned char *));
    for (i = 0; i < H; i++)
        *(a + i) = (unsigned char *) malloc(W * sizeof(unsigned char));

    ///citire matrice
    for (t = H - 1; t >= 0; t--) {
        for (j = 0; j < W; j++) {
            fread(&x, sizeof(unsigned char), 1, in);
            fread(&x, sizeof(unsigned char), 1, in);
            fread(&x, sizeof(unsigned char), 1, in);
            a[t][j] = x;
        }
        fseek(in, padd, SEEK_CUR);
    }
    fclose(in);

    ///Sablonul S
    FILE *ins = fopen(graySablon, "rb");
    if (ins == NULL) {
        printf("Nu s-a deschis sablonul S");
        exit(0);
    }

    ///preluam date
    fseek(ins, 18, SEEK_SET);
    fread(&w, sizeof(unsigned int), 1, ins);
    fread(&h, sizeof(unsigned int), 1, ins);
    fseek(ins, 54, SEEK_SET);

    ///padding sablon
    unsigned int pad = 0;
    unsigned char zero = '0';
    if ((3 * w) % 4 > 0)
        pad = 4 - (3 * w) % 4;

    ///alocarea lui S
    unsigned char **s;
    unsigned int n = w * h;///nr de pixeli
    double medf, meds = 0;

    s = (unsigned char **) malloc(h * sizeof(unsigned char *));
    for (i = 0; i < h; i++)
        *(s + i) = (unsigned char *) malloc(w * sizeof(unsigned char));

    ///citirea si calculul mediei lui S
    for (t = h - 1; t >= 0; t--) {
        for (j = 0; j < w; j++) {
            fread(&x, sizeof(unsigned char), 1, ins);
            fread(&x, sizeof(unsigned char), 1, ins);
            fread(&x, sizeof(unsigned char), 1, ins);
            s[t][j] = x;
            meds += x;
        }
        fseek(in, pad, SEEK_CUR);
    }
    meds /= n;
    fclose(ins);

    double devs = 0, devf = 0, corr = 0, val, dev;

    for (i = 0; i < h; i++)
        for (j = 0; j < w; j++)
            devs += (s[i][j] - meds) * (s[i][j] - meds);
    devs = devs / (n - 1);
    devs = sqrt(devs);

    unsigned int p, k;

    for (i = 0; i <= H - h; i++)
        for (j = 0; j <= W - w; j++) { ///pixel (i,j)
            corr = 0;
            ///media din fereastra
            medf = 0;
            for (k = i; k < h + i; k++)
                for (p = j; p < w + j; p++)
                    medf += a[k][p];
            medf /= n;

            devf = 0;
            ///deviatia din fereastra
            for (k = i; k < h + i; k++)
                for (p = j; p < w + j; p++)
                    devf += (a[k][p] - medf) * (a[k][p] - medf);
            devf = devf / (n - 1);
            devf = sqrt(devf);

            dev = devs * devf;

            for (k = i; k < h + i; k++)
                for (p = j; p < w + j; p++) {
                    val = (a[k][p] - medf) * (s[k - i][p - j] - meds);
                    val /= dev;
                    corr += val;
                }
            corr = corr / n;
            if (corr >= ps)
                adaugare(nv, v, i, j, corr, C);

        }


    for (i = 0; i < h; i++)
        free(s[i]);
    free(s);
    for (i = 0; i < H; i++)
        free(a[i]);
    free(a);
}

void pattern_recognition(char *sursa, unsigned int *n, struct detectii **v) {
    char sablon[20];
    double ps = 0.50;
    struct pixel C;

    FILE *in = fopen("fisier_tm.txt", "r");
    fscanf(in, "%s", sablon);
    fscanf(in, "%s", sablon);

    ///punem in v toate detectiile gasite
    C.b = 255;
    C.g = 0;
    C.r = 0;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 255;
    C.g = 255;
    C.r = 0;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 0;
    C.g = 255;
    C.r = 0;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 0;
    C.g = 255;
    C.r = 255;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 255;
    C.g = 0;
    C.r = 255;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 0;
    C.g = 0;
    C.r = 255;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 192;
    C.g = 192;
    C.r = 192;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 255;
    C.g = 140;
    C.r = 0;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 128;
    C.g = 0;
    C.r = 128;
    template_matching(sursa, sablon, ps, n, v, C);
    fscanf(in, "%s", sablon);
    C.b = 128;
    C.g = 0;
    C.r = 0;
    template_matching(sursa, sablon, ps, n, v, C);

    fclose(in);
}


int comp(const void *a, const void *b) {
    struct detectii va = *(struct detectii *) a;
    struct detectii vb = *(struct detectii *) b;
    if (va.corre < vb.corre) return 1;
    if (va.corre > vb.corre) return -1;
    return 0;
}


double suprapunere(struct detectii di, struct detectii dj) {
    double aria, S = 0;
    unsigned int Sdi, Sdj, h, w;
    int height, width;
    h = 15;
    w = 11;
    Sdi = Sdj = h * w;

    if (di.fi < dj.fi)
        if (di.fj < dj.fj) {
            height = (di.fi + h - dj.fi);
            width = (di.fj + w - dj.fj);
            if (width > 0 && height > 0)
                S = height * width;
            else S = 0;
        } else {
            height = (di.fi + h - dj.fi);
            width = (dj.fj + w - di.fj);
            if (width > 0 && height > 0)
                S = width * height;
            else S = 0;
        }
    else if (dj.fj < di.fj) {
        height = (dj.fi + h - di.fi);
        width = (dj.fj + w - di.fj);
        if (width > 0 && height > 0)
            S = height * width;
        else S = 0;
    } else {
        height = (dj.fi + h - di.fi);
        width = (di.fj + w - dj.fj);
        if (width > 0 && height > 0)
            S = width * height;
        else S = 0;
    }

    aria = (double) S / (Sdi + Sdj - S);
    return aria;
}


void non_maxime(unsigned int *n, struct detectii **D) {
    unsigned int i, j, na;
    struct detectii *aux;
    double x;

    for (i = 0; i < *n - 1; i++) {  //printf("%u\n",*n);
        for (j = i + 1; j < *n; j++) {
            x = suprapunere((*D)[i], (*D)[j]);
            // printf("%.3f\n",x);
            if (x > 0.2)
                //adaug(&na,&aux,(*D)[j]);
            {
                int k;
                for (k = j; k < *n - 1; k++)
                    (*D)[k] = (*D)[k + 1];
                (*n)--;
                j--;
            }
        }
    }

}

int main() {
    char caleImagine[20];
    char caleImgExtern[20];
    char keyText[20];

    ///criptarea
    FILE *in = fopen("fisier_criptare.txt", "r");
    if (in == NULL) {
        printf("Nu s-a deschis fisierul in caare se afla numele imaginilor si a cheii secrete");
        return 0;
    }
    fscanf(in, "%s%s%s", caleImagine, caleImgExtern, keyText);
    fclose(in);

    criptare(caleImagine, caleImgExtern, keyText);

    ///decriptarea
    in = fopen("fisier_decriptare.txt", "r");
    fscanf(in, "%s%s%s", caleImagine, caleImgExtern, keyText);
    fclose(in);
    decriptare(caleImagine, caleImgExtern, keyText);

    ///testul chi patrat
    printf("Valorile testului chi patrat pentru imaginea initiala sunt:\n");
    chi__patrat(caleImagine);
    printf("\n");

    printf("Valorile testului chi patrat pentru imaginea criptata sunt:\n");
    chi__patrat(caleImgExtern);
    printf("\n");

    char sursa[20];
    unsigned int n = 0;
    struct detectii *v = NULL;

    in = fopen("fisier_tm.txt", "r");
    fscanf(in, "%s", sursa);
    fclose(in);

    ///pune in v toate detectiile gasite in imaginea I
    pattern_recognition(sursa, &n, &v);

    ///sortam vectorul
    qsort(v, n, sizeof(struct detectii), comp);

    ///aplicam algoritmul de suprimare a non-maximelor
    non_maxime(&n, &v);

    ///coloram
    colorare(sursa, 15, 11, n, v);
    return 0;
}
