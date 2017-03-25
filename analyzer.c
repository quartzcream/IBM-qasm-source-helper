#include <assert.h>
#include <complex.h>
#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct complex_matrix {
  int size;
  complex double *mat;
};

void complex_mat_print(const struct complex_matrix *a, FILE *out) {
  for (int i = 0; i < a->size; ++i) {
    fprintf(out, "\n\n");
    for (int j = 0; j < a->size; ++j) {
      fprintf(out, "   %+.2f%+.2fi", creal(a->mat[i * a->size + j]),
	      cimag(a->mat[i * a->size + j]));
    }
    fprintf(out, "\n\n");
  }
}

void complex_mat_print_mag(const struct complex_matrix *a, FILE *out) {
  for (int i = 0; i < a->size; ++i) {
    fprintf(out, "\n\n");
    for (int j = 0; j < a->size; ++j) {
      fprintf(out, "   %.2f", cabs(a->mat[i * a->size + j]));
    }
    fprintf(out, "\n\n");
  }
}

void complex_mat_print_ang(const struct complex_matrix *a, FILE *out) {
  for (int i = 0; i < a->size; ++i) {
    fprintf(out, "\n\n");
    for (int j = 0; j < a->size; ++j) {
      fprintf(out, "   %.2f", carg(a->mat[i * a->size + j]));
    }
    fprintf(out, "\n\n");
  }
}

void complex_mat_to_id(struct complex_matrix *a) {
  for (int i = 0; i < a->size; ++i) {
    for (int j = 0; j < a->size; ++j) {
      a->mat[a->size * i + j] = (i == j);
    }
  }
}

// c = a*b
void complex_mat_mult(const struct complex_matrix *a,
		      const struct complex_matrix *b,
		      struct complex_matrix *c) {
  assert(a->size == b->size && a->size == c->size);
  for (int i = 0; i < a->size; ++i) {
    for (int j = 0; j < b->size; ++j) {
      c->mat[i * c->size + j] = 0;
      for (int k = 0; k < c->size; ++k) {
	c->mat[i * c->size + j] +=
	    a->mat[i * a->size + k] * b->mat[k * b->size + j];
      }
    }
  }
  return;
}

struct named_gate {
  char *name;
  struct complex_matrix mat;
};

int avail_gate_c;
struct named_gate *avail_gates;

void init_avail_gates() {
  avail_gates = malloc(50 * sizeof(struct named_gate));
  int i = 0;
  {
    avail_gates[i].name = "qreg";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {1, 0, 0, 1};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "id";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {1, 0, 0, 1};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "x";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {0, 1, 1, 0};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "y";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {0, -I, I, 0};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "z";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {1, 0, 0, -1};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "h";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {sqrt(0.5), sqrt(0.5), sqrt(0.5), -sqrt(0.5)};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "s";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {1, 0, 0, I};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "sdg";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {1, 0, 0, -I};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "t";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {1, 0, 0, (1 + I) / sqrt(2.0)};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "tdg";
    avail_gates[i].mat.size = 2;
    avail_gates[i].mat.mat = malloc(2 * 2 * sizeof(complex double));
    complex double tmp[4] = {1, 0, 0, (1 - I) / sqrt(2.0)};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  {
    avail_gates[i].name = "cx";
    avail_gates[i].mat.size = 4;
    avail_gates[i].mat.mat = malloc(4 * 4 * sizeof(complex double));
    complex double tmp[16] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0};
    memcpy(avail_gates[i].mat.mat, tmp, sizeof(tmp));
    ++i;
  }
  avail_gate_c = i;
}

void clean_avail_gates() {
  for (int i = 0; i < avail_gate_c; ++i) {
    free(avail_gates[i].mat.mat);
  }
  free(avail_gates);
}

struct gate_pos {
  int len;
  int id;
  int *ord;
};

int read_gate(const char *gate_name, char *ord_str, struct gate_pos *gate_p) {
  gate_p->len = 0;
  for (int i = 0; i < avail_gate_c; ++i) {
    if (!strcmp(gate_name, avail_gates[i].name)) {
      gate_p->id = i;
      gate_p->len = __builtin_ctz(avail_gates[i].mat.size);
      gate_p->ord = malloc(gate_p->len * sizeof(int));
      break;
    }
  }
  if (gate_p->len == 0) {
    return 1;
  }
  char *token = strtok(ord_str, "[]");
  for (int i = 0; i < gate_p->len;) {
    assert(token);
    char *tmp;
    int val = strtol(token, &tmp, 10);
    if (*tmp == '\0') {
      gate_p->ord[i] = val;
      ++i;
    }
    token = strtok(NULL, "[]");
  }
  return 0;
}

void gate_to_trans_mat(const struct gate_pos *gate,
		       struct complex_matrix *mat) {
  assert(__builtin_popcount(mat->size) == 1);
  for (int i = 0; i < mat->size; ++i) {
    for (int j = 0; j < mat->size; ++j) {
      int gate_i = 0, gate_j = 0;
      int rest_i = i, rest_j = j;
      for (int part = 0; part < gate->len; ++part) {
	gate_i |= ((i & (1 << gate->ord[part])) && 1) << part;
	rest_i &= ~(1 << gate->ord[part]);

	gate_j |= ((j & (1 << gate->ord[part])) && 1) << part;
	rest_j &= ~(1 << gate->ord[part]);
      }
      if (rest_i == rest_j) {
	mat->mat[i * mat->size + j] =
	    avail_gates[gate->id]
		.mat.mat[gate_i * avail_gates[gate->id].mat.size + gate_j];
      } else {
	mat->mat[i * mat->size + j] = 0;
      }
    }
  }
}

int main(int argc, char **argv) {
  init_avail_gates();

  char *file_name = "inp.qasm";
  if (argc > 1) {
    file_name = argv[1];
  }
  FILE *f_source = fopen(file_name, "r");
  fseek(f_source, 0, SEEK_END);
  long len = ftell(f_source);
  fseek(f_source, 0, SEEK_SET);
  char *source = malloc(len + 1);
  fread(source, 1, len, f_source);
  source[len] = '\0';
  fclose(f_source);

  size_t line_l = 10;
  size_t line_c = 0;
  char **lines = malloc(line_l * sizeof(char *));
  char *token = strtok(source, "\n\r");
  for (; token; ++line_c) {
    if (line_c == line_l) {
      lines = realloc(lines, line_l * 2 * sizeof(char *));
      line_l *= 2;
    }
    lines[line_c] = token;
    token = strtok(NULL, "\n\r");
  }

  char *source_bak = malloc(len + 1);
  memcpy(source_bak, source, len + 1);

  int qubit_c = 0;
  struct gate_pos *gate_ord = malloc(line_c * sizeof(struct gate_pos));
  for (int i = 0; i < line_c; ++i) {
    char *gate_name = strtok(lines[i], " \t\n\v\f\r");
    char *ord_str = strtok(NULL, " \t\n\v\f\r");

    int empty_line = read_gate(gate_name, ord_str, gate_ord + i);
    if (empty_line) {
      gate_ord[i].id = 1;
      gate_ord[i].len = 1;
      gate_ord[i].ord = calloc(1, sizeof(int));
    }
    if (gate_ord[i].id == 0) {
      assert(qubit_c == 0);
      qubit_c = gate_ord[i].ord[0];
      gate_ord[i].id = 1;
      gate_ord[i].len = 1;
      gate_ord[i].ord[0] = 0;
    }
  }

  assert(qubit_c != 0);
  struct complex_matrix *trans_ord =
      malloc(line_c * sizeof(struct complex_matrix));
  for (int i = 0; i < line_c; ++i) {
    trans_ord[i].size = 1 << qubit_c;
    trans_ord[i].mat =
	malloc((1 << qubit_c) * (1 << qubit_c) * sizeof(complex double));
    gate_to_trans_mat(gate_ord + i, trans_ord + i);
  }

  /*
	for(int i = 0; i < line_c; ++i){
		printf("%d:\n", i);
		complex_mat_print(trans_ord + i, stdout);
	}
*/

  char buffer[100];
  struct complex_matrix mat_out = {
      1 << qubit_c,
      malloc((1 << qubit_c) * (1 << qubit_c) * sizeof(complex double))};
  struct complex_matrix tmp = {
      1 << qubit_c,
      malloc((1 << qubit_c) * (1 << qubit_c) * sizeof(complex double))};
  while (scanf("%s", buffer) == 1) {
    printf(">>>");
    fflush(stdout);
    if (!strcmp(buffer, "p") || !strcmp(buffer, "abs") ||
	!strcmp(buffer, "ang")) {
      int beg, end;
      scanf("%d %d", &beg, &end);
      if (beg < 0 || beg > end || end > line_c) {
	break;
      }
      complex_mat_to_id(&mat_out);
      for (int i = beg; i < end; ++i) {
	complex_mat_mult(trans_ord + i, &mat_out, &tmp);
	complex double *p_tmp = tmp.mat;
	tmp.mat = mat_out.mat;
	mat_out.mat = p_tmp;
      }
      printf("\n");
      for (int i = beg; i < end; ++i) {
	printf("%s\n", source_bak + (lines[i] - source));
      }
      if (!strcmp(buffer, "p")) {
	complex_mat_print(&mat_out, stdout);
      } else if (!strcmp(buffer, "abs")) {
	complex_mat_print_mag(&mat_out, stdout);
      } else if (!strcmp(buffer, "ang")) {
	complex_mat_print_ang(&mat_out, stdout);
      }
      continue;
    }
    printf("\n");
    break;
  }

  free(mat_out.mat);
  free(tmp.mat);

  for (int i = 0; i < line_c; ++i) {
    free(gate_ord[i].ord);
  }
  free(gate_ord);

  for (int i = 0; i < line_c; ++i) {
    free(trans_ord[i].mat);
  }
  free(trans_ord);

  free(lines);
  free(source_bak);
  free(source);
  clean_avail_gates();
}
