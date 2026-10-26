# ==============================================================
# TRABALHO PRÁTICO II — TESTES DE HIPÓTESES PARAMÉTRICOS
# ==============================================================

library(tidyverse)
library(car)

# ==============================================================
# 1. IMPORTAÇÃO E PREPARAÇÃO
# ==============================================================

# Se Dados já estiver carregado, comenta a linha seguinte.
Dados <- read.delim(file.choose(), check.names = FALSE)

# Confirmar se as colunas necessárias existem
colunas_necessarias <- c(
  "ID", "TREAT", "AGE", "WEIGHT",
  "LIVINGSTATUS", "SMOKESTATUS", "GENDER", "VISIT"
)

faltam <- setdiff(colunas_necessarias, names(Dados))

if (length(faltam) > 0) {
  stop(paste("Faltam estas colunas:", paste(faltam, collapse = ", ")))
}

# Limpar texto e garantir que idade e peso são numéricos
Dados <- Dados |>
  mutate(
    across(
      all_of(c("TREAT", "LIVINGSTATUS", "SMOKESTATUS", "GENDER", "VISIT")),
      ~ toupper(trimws(as.character(.x)))
    ),
    AGE = as.numeric(AGE),
    WEIGHT = as.numeric(WEIGHT)
  )

# Cada participante deve ter apenas uma observação PRE e uma POST
if (anyDuplicated(Dados[c("ID", "VISIT")]) > 0) {
  stop("Há mais de uma observação para o mesmo participante e visita.")
}

# Transformar para uma linha por participante
Dados_largo <- Dados |>
  pivot_wider(
    id_cols = c(ID, TREAT, AGE, LIVINGSTATUS, SMOKESTATUS, GENDER),
    names_from = VISIT,
    values_from = WEIGHT
  )

if (!all(c("PRE", "POST") %in% names(Dados_largo))) {
  stop("Não foram encontradas as visitas PRE e POST.")
}

if (anyNA(Dados_largo$PRE) || anyNA(Dados_largo$POST)) {
  stop("Existem participantes sem peso PRE ou POST.")
}

# Criar tercis etários
limites_idade <- as.numeric(
  quantile(
    Dados_largo$AGE,
    probs = c(0, 1/3, 2/3, 1),
    na.rm = TRUE
  )
)

if (length(unique(limites_idade)) != 4) {
  stop("Os limites dos tercis etários repetem-se; é necessário rever as faixas.")
}

# ==============================================================
# 2. CRIAÇÃO E RECODIFICAÇÃO DAS VARIÁVEIS
# ==============================================================

Dados_largo <- Dados_largo |>
  mutate(
    # Valor positivo significa perda de peso
    PERDA = PRE - POST,
    PERDA_PERC = 100 * PERDA / PRE,
    PERDA_5 = PERDA_PERC >= 5,
    
    # Manter os três braços para análises complementares
    GRUPO3 = factor(
      case_when(
        TREAT == "CONTROL"   ~ "Controlo",
        TREAT == "PLACEBO"   ~ "Placebo",
        TREAT == "TREATMENT" ~ "Tratamento",
        TRUE ~ NA_character_
      ),
      levels = c("Controlo", "Placebo", "Tratamento")
    ),
    
    # Categoria binária pedida:
    # Controlo = CONTROL + PLACEBO
    GRUPO2 = factor(
      case_when(
        TREAT %in% c("CONTROL", "PLACEBO") ~ "Controlo",
        TREAT == "TREATMENT"               ~ "Tratamento",
        TRUE ~ NA_character_
      ),
      levels = c("Controlo", "Tratamento")
    ),
    
    GENERO = factor(
      case_when(
        GENDER == "FEMALE" ~ "Feminino",
        GENDER == "MALE"   ~ "Masculino",
        TRUE ~ NA_character_
      ),
      levels = c("Feminino", "Masculino")
    ),
    
    # Fumador actual versus não fumador actual.
    # QUIT-SMOKER é classificado como não fumador actual.
    FUMADOR2 = factor(
      case_when(
        SMOKESTATUS == "SMOKER" ~ "Fumador",
        SMOKESTATUS %in% c("NON-SMOKER", "QUIT-SMOKER") ~ "Não fumador",
        TRUE ~ NA_character_
      ),
      levels = c("Não fumador", "Fumador")
    ),
    
    # Esta variável representa coabitação, não estado civil.
    VIVE2 = factor(
      case_when(
        str_detect(LIVINGSTATUS, "^ALONE") ~ "Vive sozinho(a)",
        str_detect(LIVINGSTATUS, "^WITH")  ~ "Vive acompanhado(a)",
        TRUE ~ NA_character_
      ),
      levels = c("Vive acompanhado(a)", "Vive sozinho(a)")
    ),
    
    # Faixas etárias baseadas em tercis
    FAIXA_IDADE = cut(
      AGE,
      breaks = limites_idade,
      include.lowest = TRUE,
      labels = c("Tercil 1", "Tercil 2", "Tercil 3"),
      ordered_result = TRUE
    ),
    
    # Indicador opcional de coabitação com parceiro(a).
 
    PARCEIRO2 = factor(
      case_when(
        LIVINGSTATUS %in% c(
          "WITH PARTNER",
          "WITH PARTNER AND CHILDREN"
        ) ~ "Coabita com parceiro(a)",
        
        LIVINGSTATUS %in% c(
          "ALONE (SINGLE)",
          "ALONE (DIVORCED)",
          "ALONE (WIDOWER)",
          "WITH FRIENDS"
        ) ~ "Sem parceiro(a) indicado",
        
        TRUE ~ NA_character_
      ),
      levels = c("Sem parceiro(a) indicado", "Coabita com parceiro(a)")
    )
  )

# Confirmar as recodificações
table(Dados_largo$TREAT, Dados_largo$GRUPO3, useNA = "ifany")
table(Dados_largo$TREAT, Dados_largo$GRUPO2, useNA = "ifany")
table(Dados_largo$SMOKESTATUS, Dados_largo$FUMADOR2, useNA = "ifany")
table(Dados_largo$LIVINGSTATUS, Dados_largo$VIVE2, useNA = "ifany")

# No conjunto de dados apresentado, espera-se:
# GRUPO2: 80 Controlo (CONTROL + PLACEBO) e 40 Tratamento
# FUMADOR2: 61 Não fumadores e 59 Fumadores actuais
# VIVE2: 62 acompanhados e 58 sozinhos
table(Dados_largo$GRUPO2)
table(Dados_largo$FUMADOR2)
table(Dados_largo$VIVE2)

# Verificar se há categorias não recodificadas
variaveis_verificar <- c(
  "GRUPO3", "GRUPO2", "GENERO",
  "FUMADOR2", "VIVE2", "FAIXA_IDADE"
)

if (anyNA(Dados_largo[variaveis_verificar])) {
  stop("Há valores em falta numa ou mais variáveis recodificadas.")
}

# ==============================================================
# 3. FUNÇÃO PARA O TAMANHO DE EFEITO
# ==============================================================

# d de Cohen absoluto, apenas como medida complementar
cohen_d_abs <- function(x, g) {
  completos <- complete.cases(x, g)
  x <- x[completos]
  g <- droplevels(factor(g[completos]))
  
  if (nlevels(g) != 2) {
    return(NA_real_)
  }
  
  nivel1 <- levels(g)[1]
  nivel2 <- levels(g)[2]
  
  x1 <- x[g == nivel1]
  x2 <- x[g == nivel2]
  
  if (length(x1) < 2 || length(x2) < 2) {
    return(NA_real_)
  }
  
  dp_agregado <- sqrt(
    (
      (length(x1) - 1) * stats::var(x1) +
        (length(x2) - 1) * stats::var(x2)
    ) /
      (length(x1) + length(x2) - 2)
  )
  
  if (dp_agregado == 0) {
    return(NA_real_)
  }
  
  abs(mean(x1) - mean(x2)) / dp_agregado
}

# ==============================================================
# 4. ANÁLISE DESCRITIVA
# ==============================================================

nrow(Dados_largo)

# Resumo dos três braços originais
resumo_3_grupos <- Dados_largo |>
  group_by(GRUPO3) |>
  summarise(
    n = n(),
    idade_media = mean(AGE),
    peso_PRE_medio = mean(PRE),
    peso_POST_medio = mean(POST),
    perda_media = mean(PERDA),
    perda_mediana = median(PERDA),
    perda_dp = sd(PERDA),
    .groups = "drop"
  )

print(resumo_3_grupos)

# Resumo da comparação binária pedida
resumo_2_grupos <- Dados_largo |>
  group_by(GRUPO2) |>
  summarise(
    n = n(),
    peso_PRE_medio = mean(PRE),
    peso_POST_medio = mean(POST),
    perda_media = mean(PERDA),
    perda_mediana = median(PERDA),
    perda_dp = sd(PERDA),
    .groups = "drop"
  )

print(resumo_2_grupos)

# Gráfico da perda por braço original
ggplot(Dados_largo, aes(x = GRUPO3, y = PERDA, fill = GRUPO3)) +
  geom_boxplot() +
  labs(
    title = "Perda de peso por grupo",
    x = "Grupo",
    y = "Perda de peso (kg)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# ==============================================================
# 5. VERIFICAÇÃO DOS PRESSUPOSTOS
# ==============================================================

# Normalidade da perda dentro de cada grupo
normalidade_grupos <- function(dados, grupo) {
  dados |>
    group_by(.data[[grupo]]) |>
    summarise(
      n = n(),
      W_shapiro = unname(shapiro.test(PERDA)$statistic),
      p_shapiro = shapiro.test(PERDA)$p.value,
      .groups = "drop"
    )
}

cat("\nNormalidade nos três grupos:\n")
print(normalidade_grupos(Dados_largo, "GRUPO3"))

cat("\nNormalidade nos dois grupos:\n")
print(normalidade_grupos(Dados_largo, "GRUPO2"))

# Gráficos Q-Q da perda por grupo
ggplot(Dados_largo, aes(sample = PERDA)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ GRUPO3) +
  labs(
    title = "Gráficos Q-Q da perda de peso",
    x = "Quantis teóricos",
    y = "Quantis observados"
  ) +
  theme_minimal()

# Homogeneidade das variâncias
cat("\nLevene para os três grupos:\n")
print(car::leveneTest(PERDA ~ GRUPO3, data = Dados_largo, center = median))

cat("\nLevene para os dois grupos:\n")
print(car::leveneTest(PERDA ~ GRUPO2, data = Dados_largo, center = median))

# ==============================================================
# 6. OBJECTIVO 1 — EFEITO DO TRATAMENTO
# ==============================================================

# 6.1 Teste t emparelhado: PRE versus POST em cada braço original
teste_pareado_por_grupo <- function(dados, nome_grupo) {
  dados |>
    group_by(.data[[nome_grupo]]) |>
    group_modify(~ {
      d <- .x
      tt <- t.test(d$PRE, d$POST, paired = TRUE)
      
      tibble(
        n = nrow(d),
        perda_media = mean(d$PERDA),
        perda_dp = sd(d$PERDA),
        t = unname(tt$statistic),
        gl = unname(tt$parameter),
        valor_p = tt$p.value,
        IC95_inf = tt$conf.int[1],
        IC95_sup = tt$conf.int[2],
        d_z = mean(d$PERDA) / sd(d$PERDA)
      )
    }) |>
    ungroup()
}

cat("\nTestes t emparelhados por braço:\n")
resultado_pareado <- teste_pareado_por_grupo(Dados_largo, "GRUPO3")
print(resultado_pareado)

# Nota: estes testes avaliam PRE versus POST dentro de cada braço.
# A comparação entre os braços é feita abaixo, usando PERDA.

# 6.2 Comparação principal pedida: Controlo (CONTROL + PLACEBO) vs Tratamento
# t.test usa, por padrão, o teste de Welch (var.equal = FALSE).
teste_grupo2 <- t.test(
  PERDA ~ GRUPO2,
  data = Dados_largo,
  var.equal = FALSE
)

cat("\nTeste t de Welch: Controlo (CONTROL + PLACEBO) vs Tratamento:\n")
print(teste_grupo2)

niveis_grupo2 <- levels(Dados_largo$GRUPO2)

resultado_grupo2 <- tibble(
  comparacao = paste(niveis_grupo2[1], "-", niveis_grupo2[2]),
  n1 = sum(Dados_largo$GRUPO2 == niveis_grupo2[1]),
  n2 = sum(Dados_largo$GRUPO2 == niveis_grupo2[2]),
  media1 = unname(teste_grupo2$estimate[1]),
  media2 = unname(teste_grupo2$estimate[2]),
  diferenca_media = unname(teste_grupo2$estimate[1] - teste_grupo2$estimate[2]),
  IC95_inf = teste_grupo2$conf.int[1],
  IC95_sup = teste_grupo2$conf.int[2],
  t = unname(teste_grupo2$statistic),
  gl = unname(teste_grupo2$parameter),
  valor_p = teste_grupo2$p.value,
  d_cohen_abs = cohen_d_abs(Dados_largo$PERDA, Dados_largo$GRUPO2)
)

print(resultado_grupo2)

# 6.3 Análise complementar: manter Controlo, Placebo e Tratamento separados
# Welch ANOVA, adequada quando as variâncias são diferentes
cat("\nANOVA de Welch para os três grupos:\n")
print(oneway.test(PERDA ~ GRUPO3, data = Dados_largo, var.equal = FALSE))

# Comparações duas a duas de Welch, com correcção de Holm
cat("\nComparações duas a duas de Welch, com correcção de Holm:\n")
print(
  pairwise.t.test(
    Dados_largo$PERDA,
    Dados_largo$GRUPO3,
    pool.sd = FALSE,
    p.adjust.method = "holm"
  )
)

# Comparação directa Controlo versus Tratamento, excluindo o Placebo
Dados_CT <- Dados_largo |>
  filter(GRUPO3 %in% c("Controlo", "Tratamento")) |>
  mutate(GRUPO3 = droplevels(GRUPO3))

cat("\nTeste t de Welch: apenas Controlo versus Tratamento:\n")
print(t.test(PERDA ~ GRUPO3, data = Dados_CT, var.equal = FALSE))

cat("\nd de Cohen absoluto: apenas Controlo versus Tratamento:\n")
print(cohen_d_abs(Dados_CT$PERDA, Dados_CT$GRUPO3))

# ANOVA clássica de três grupos, apenas como análise complementar.
# Se Levene for significativo, a referência deve ser a ANOVA de Welch acima.
modelo_classico_3 <- aov(PERDA ~ GRUPO3, data = Dados_largo)

cat("\nANOVA clássica de três grupos (análise complementar):\n")
print(summary(modelo_classico_3))

somas_quadrados <- summary(modelo_classico_3)[[1]][["Sum Sq"]]
eta2_classico <- somas_quadrados[1] / sum(somas_quadrados)

cat("\nEta-quadrado da ANOVA clássica:\n")
print(eta2_classico)

# ==============================================================
# 7. OBJECTIVO 2 — FACTORES ASSOCIADOS À PERDA NO TRATAMENTO
# ==============================================================

tratamento <- Dados_largo |>
  filter(GRUPO3 == "Tratamento")

# Analisar estes factores dentro do grupo Tratamento
variaveis_2_niveis <- c("GENERO", "FUMADOR2", "VIVE2")

# Se o docente aceitar coabitação com parceiro(a) como aproximação
# ao estado marital, pode substituir VIVE2 por PARCEIRO2:
# variaveis_2_niveis <- c("GENERO", "FUMADOR2", "PARCEIRO2")

# Descritivas e normalidade por categoria
descrever_factor <- function(dados, nome_var) {
  tab <- dados |>
    group_by(.data[[nome_var]]) |>
    summarise(
      n = n(),
      perda_media = mean(PERDA),
      perda_dp = sd(PERDA),
      W_shapiro = unname(shapiro.test(PERDA)$statistic),
      p_shapiro = shapiro.test(PERDA)$p.value,
      .groups = "drop"
    )
  
  names(tab)[1] <- "categoria"
  
  tab |>
    mutate(variavel = nome_var, .before = 1)
}

cat("\nDescritivas e normalidade por categoria, no grupo Tratamento:\n")
descritivas_factores <- map_dfr(
  variaveis_2_niveis,
  ~ descrever_factor(tratamento, .x)
)

print(descritivas_factores, width = Inf)

# Função para extrair o resultado do teste de Levene
extrair_levene <- function(dados, nome_grupo) {
  formula <- reformulate(nome_grupo, response = "PERDA")
  lev <- car::leveneTest(formula, data = dados, center = median)
  
  tibble(
    variavel = nome_grupo,
    gl1 = as.numeric(lev[1, "Df"]),
    gl2 = as.numeric(lev[2, "Df"]),
    F = as.numeric(lev[1, "F value"]),
    valor_p = as.numeric(lev[1, "Pr(>F)"])
  )
}

cat("\nTestes de Levene para os factores no grupo Tratamento:\n")
levene_factores <- map_dfr(
  variaveis_2_niveis,
  ~ extrair_levene(tratamento, .x)
)

print(levene_factores)

# Função para o teste t de Welch entre duas categorias
teste_welch_factor <- function(dados, nome_var) {
  d <- dados |>
    select(PERDA, all_of(nome_var)) |>
    drop_na()
  
  d[[nome_var]] <- droplevels(factor(d[[nome_var]]))
  niveis <- levels(d[[nome_var]])
  
  if (length(niveis) != 2) {
    stop(paste("A variável", nome_var, "não tem exactamente duas categorias."))
  }
  
  formula <- reformulate(nome_var, response = "PERDA")
  tt <- t.test(formula, data = d, var.equal = FALSE)
  
  tibble(
    variavel = nome_var,
    categoria1 = niveis[1],
    categoria2 = niveis[2],
    n1 = sum(d[[nome_var]] == niveis[1]),
    n2 = sum(d[[nome_var]] == niveis[2]),
    media1 = unname(tt$estimate[1]),
    media2 = unname(tt$estimate[2]),
    diferenca_media_1_menos_2 =
      unname(tt$estimate[1] - tt$estimate[2]),
    IC95_inf = tt$conf.int[1],
    IC95_sup = tt$conf.int[2],
    t = unname(tt$statistic),
    gl = unname(tt$parameter),
    valor_p = tt$p.value,
    d_cohen_abs = cohen_d_abs(d$PERDA, d[[nome_var]])
  )
}

cat("\nTestes t de Welch dos factores no grupo Tratamento:\n")
resultado_factores <- map_dfr(
  variaveis_2_niveis,
  ~ teste_welch_factor(tratamento, .x)
) |>
  mutate(
    p_Holm = p.adjust(valor_p, method = "holm")
  )

print(resultado_factores, width = Inf)

# ==============================================================
# 8. IDADE NO GRUPO TRATAMENTO
# ==============================================================

cat("\nPerda de peso por faixa etária no grupo Tratamento:\n")
print(
  tratamento |>
    group_by(FAIXA_IDADE) |>
    summarise(
      n = n(),
      perda_media = mean(PERDA),
      perda_dp = sd(PERDA),
      .groups = "drop"
    )
)

cat("\nLevene para as faixas etárias no grupo Tratamento:\n")
levene_idade <- car::leveneTest(
  PERDA ~ FAIXA_IDADE,
  data = tratamento,
  center = median
)
print(levene_idade)

# ANOVA clássica de um factor
modelo_idade <- aov(PERDA ~ FAIXA_IDADE, data = tratamento)

cat("\nNormalidade dos resíduos da ANOVA da idade:\n")
print(shapiro.test(residuals(modelo_idade)))

cat("\nANOVA clássica da idade:\n")
print(summary(modelo_idade))

# ANOVA de Welch e comparações de Welch como análise complementar
cat("\nANOVA de Welch para as faixas etárias:\n")
print(oneway.test(PERDA ~ FAIXA_IDADE, data = tratamento, var.equal = FALSE))

cat("\nComparações de Welch entre faixas etárias, com correcção de Holm:\n")
print(
  pairwise.t.test(
    tratamento$PERDA,
    tratamento$FAIXA_IDADE,
    pool.sd = FALSE,
    p.adjust.method = "holm"
  )
)

# Tukey só deve ser interpretado se os pressupostos da ANOVA clássica
# forem razoáveis e a ANOVA global for significativa.

