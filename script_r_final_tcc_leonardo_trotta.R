# TCC — MBA Data Science e Analytics (USP/ESALQ)
# Consumo de Suplementos Físicos de RPG no Brasil
# Script de análise — Leonardo Trotta Banci
# Entrada: dados_tratados.csv

library(dplyr)
library(readr)
library(psych)
library(GPArotation)
library(car)
library(MASS)
library(ResourceSelection)
library(pROC)
library(performance)
library(caret)

###############################################################################
# ETAPA 1 — Base tratada e estatísticas descritivas
###############################################################################

# 1.1 Importação e verificações
arquivo_base <- "dados_tratados.csv"

df <- read_csv2(arquivo_base, locale = locale(encoding = "UTF-8")) #mude para read_csv de acordo com a configuração do seu dispositivo

# 1.2 Tipagem das variáveis
df <- df %>%
  mutate(
    genero = factor(genero),
    
    renda = factor(renda, levels = c(
      "Até 1 salário mínimo", "De 1 a 3 salários mínimos",
      "De 3 a 5 salários mínimos", "De 5 a 10 salários mínimos",
      "Acima de 10 salários mínimos"
    ), ordered = TRUE),
    
    escolaridade = factor(escolaridade, levels = c(
      "Ensino médio incompleto", "Ensino médio completo",
      "Ensino superior incompleto", "Ensino superior completo",
      "Pós-graduação (especialização, mestrado ou doutorado)"
    ), ordered = TRUE),
    
    tempo_rpg = factor(tempo_rpg, levels = c(
      "Até 1 ano", "Entre 1 e 2 anos", "Entre 2 e 3 anos",
      "Entre 3 e 4 anos", "5 anos ou mais"
    ), ordered = TRUE),
    
    freq_sessoes = factor(freq_sessoes, levels = c(
      "Não estou jogando no momento", "Raramente", "Mensal",
      "Quinzenal", "Semanal"
    ), ordered = TRUE),
    
    compra_pdf = factor(compra_pdf, levels = c(
      "Não utilizo PDFs", "Sim, mas apenas gratuitos/pirataria",
      "Sim, compro regularmente"
    )),
    
    papel_mesa             = factor(papel_mesa),
    sistema_rpg            = factor(sistema_rpg),
    consome_youtube_twitch = factor(consome_youtube_twitch),
    zona_residencia        = factor(zona_residencia),
    consumidor             = factor(consumidor, levels = c("Não comprador", "Comprador")),
    
    across(c(TAM_1, TAM_2, TAM_3, EST_1, EST_2, EST_3,
             SUBNORM_1, AUT_1, SUBNORM_2, VAM_1, VAM_2, VAM_3, SUBNORM_3),
           ~ factor(.x, levels = c(
             "Discordo totalmente", "Discordo", "Discordo parcialmente",
             "Neutro", "Concordo parcialmente", "Concordo", "Concordo totalmente"
           ), ordered = TRUE))
  )

# 1.3 Verificação de valores ausentes após a tipagem
usadas <- c("consumidor", "renda", "freq_sessoes", "papel_mesa", "compra_pdf",
            "TAM_1", "TAM_3", "EST_1", "EST_3", "SUBNORM_2", "VAM_2", "VAM_3")
na_usadas <- colSums(is.na(df[, usadas]))
if (any(na_usadas > 0)) {
  stop("Valores ausentes após a tipagem em: ",
       paste(names(na_usadas)[na_usadas > 0], collapse = ", "),
       ". Confira se os rótulos do CSV batem com os níveis declarados.")
}
na_outras <- colSums(is.na(df[, c("genero", "escolaridade", "tempo_rpg", "sistema_rpg",
                                  "consome_youtube_twitch", "zona_residencia")]))
if (any(na_outras > 0)) {
  warning("Valores ausentes após a tipagem em: ",
          paste(names(na_outras)[na_outras > 0], collapse = ", "))
}

# 1.4 Estatísticas descritivas
pct <- function(x) round(100 * prop.table(table(x)), 1)

pct(df$genero)
c(media = mean(df$idade, na.rm = TRUE), dp = sd(df$idade, na.rm = TRUE),
  n = sum(!is.na(df$idade)))
pct(df$escolaridade)
head(sort(pct(df$estado), decreasing = TRUE), 3)
pct(df$renda)
pct(df$papel_mesa)
pct(df$tempo_rpg)
pct(df$freq_sessoes)
pct(df$sistema_rpg)
pct(df$consome_youtube_twitch)
pct(df$compra_pdf)

pct(df$consumidor)
pct(df$qtd_suplementos_12m[df$consumidor == "Comprador"])

# 1.5 Taxa de compra por categoria
round(100 * prop.table(table(df$papel_mesa, df$consumidor), 1)[, "Comprador"], 1)
round(100 * prop.table(table(df$renda, df$consumidor), 1)[, "Comprador"], 1)
round(100 * prop.table(table(df$compra_pdf, df$consumidor), 1)[, "Comprador"], 1)

###############################################################################
# ETAPA 2 — Análise fatorial exploratória e escores dos construtos
###############################################################################

# 2.1 Itens da análise (VAM_2 e VAM_3 com escala invertida)
itens <- df %>%
  dplyr::select(TAM_1, TAM_3, EST_1, EST_3, SUBNORM_2, VAM_2, VAM_3) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(VAM_2 = 8 - VAM_2, VAM_3 = 8 - VAM_3)

blocos <- list(
  valorizacao_posse   = c("TAM_1", "TAM_3", "EST_1", "EST_3", "SUBNORM_2"),
  resistencia_digital = c("VAM_2", "VAM_3")
)

# 2.2 Adequação da amostra
poly_cor <- polychoric(itens)$rho

KMO(poly_cor)
cortest.bartlett(poly_cor, n = nrow(itens))

# 2.3 Extração e rotação
efa <- fa(poly_cor, nfactors = 2, n.obs = nrow(itens), fm = "pa", rotate = "oblimin")

print(efa$loadings, cutoff = 0.30, sort = TRUE)
print(efa$loadings, cutoff = 0)
fa.parallel(poly_cor, n.obs = nrow(itens), fa = "fa", fm = "pa", n.iter = 100)
efa$Vaccounted
efa$communality
efa$Phi

# 2.4 Confiabilidade
alpha_valorizacao <- psych::alpha(itens[, blocos$valorizacao_posse])
alpha_valorizacao$total$raw_alpha

r_vam <- cor(itens$VAM_2, itens$VAM_3)
spearman_brown <- (2 * r_vam) / (1 + r_vam)
r_vam
spearman_brown

r_poli <- poly_cor["VAM_2", "VAM_3"]
spearman_brown_poli <- (2 * r_poli) / (1 + r_poli)
c(r_pearson = r_vam, r_policorica = r_poli,
  SB_pearson = spearman_brown, SB_policorica = spearman_brown_poli)

# 2.5 Escores dos construtos
df_final <- df %>%
  mutate(
    escala_valorizacao_posse   = rowMeans(itens[, blocos$valorizacao_posse]),
    escala_resistencia_digital = rowMeans(itens[, c("VAM_2", "VAM_3")])
  )

cor(df_final$escala_valorizacao_posse, df_final$escala_resistencia_digital)

# 2.6 Médias por grupo e efeito univariado
df_final %>%
  group_by(consumidor) %>%
  summarise(across(c(escala_valorizacao_posse, escala_resistencia_digital),
                   list(media = mean, dp = sd)), .groups = "drop")

for (col in c("escala_valorizacao_posse", "escala_resistencia_digital")) {
  m <- glm(df_final$consumidor ~ df_final[[col]], family = binomial)
  cat("\n", col, "-> OR:", round(exp(coef(m))[2], 3),
      " p-valor:", format.pval(summary(m)$coefficients[2, 4], digits = 3, eps = 1e-4), "\n")
}

# 2.7 Etapa exploratória: estrutura teórica de 12 itens e cenários comparados
itens12 <- df %>%
  dplyr::select(TAM_1, TAM_2, TAM_3, EST_1, EST_2, EST_3,
                SUBNORM_1, SUBNORM_2, SUBNORM_3, VAM_1, VAM_2, VAM_3) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(VAM_2 = 8 - VAM_2, VAM_3 = 8 - VAM_3)

blocos12 <- list(
  TAM     = c("TAM_1", "TAM_2", "TAM_3"),
  EST     = c("EST_1", "EST_2", "EST_3"),
  SUBNORM = c("SUBNORM_1", "SUBNORM_2", "SUBNORM_3"),
  VAM     = c("VAM_1", "VAM_2", "VAM_3")
)

round(sapply(blocos12, function(v) psych::alpha(itens12[, v])$total$raw_alpha), 2)
round(psych::alpha(itens12[, blocos12$VAM])$alpha.drop[, "raw_alpha"], 3)

variancia_afe <- function(vars, k) {
  r <- polychoric(itens12[, vars])$rho
  f <- fa(r, nfactors = k, n.obs = nrow(itens12), fm = "pa", rotate = "oblimin")
  round(sum(f$Vaccounted["Proportion Var", ]), 3)
}

vars_cenario1 <- c("TAM_1", "TAM_2", "TAM_3", "EST_1", "EST_2", "EST_3",
                   "VAM_1", "VAM_2", "VAM_3")
vars_cenario2 <- c("TAM_1", "TAM_2", "TAM_3", "EST_1", "EST_2", "EST_3",
                   "SUBNORM_2", "VAM_1", "VAM_2", "VAM_3")
c(estrutura_12itens_4fatores = variancia_afe(names(itens12), 4),
  cenario1_3fatores          = variancia_afe(vars_cenario1, 3),
  cenario2_3fatores          = variancia_afe(vars_cenario2, 3),
  cenario2_4fatores          = variancia_afe(vars_cenario2, 4))

# 2.8 Critério de retenção do número de fatores
round(eigen(poly_cor)$values, 2)
efa3 <- fa(poly_cor, nfactors = 3, n.obs = nrow(itens), fm = "pa", rotate = "oblimin")
print(efa3$loadings, cutoff = 0)
efa3$Phi

###############################################################################
# ETAPA 3 — Regressão logística binária
###############################################################################

df <- df_final

# 3.1 Eventos por variável
eventos <- min(table(df$consumidor))
eventos

# 3.2 Codificação dos preditores
df <- df %>%
  mutate(
    renda_num = case_when(
      renda == "Até 1 salário mínimo"          ~ 1,
      renda == "De 1 a 3 salários mínimos"     ~ 2,
      renda == "De 3 a 5 salários mínimos"     ~ 3,
      renda == "De 5 a 10 salários mínimos"    ~ 4,
      renda == "Acima de 10 salários mínimos"  ~ 5
    ),
    freq_sessoes_num = case_when(
      freq_sessoes == "Não estou jogando no momento" ~ 1,
      freq_sessoes == "Raramente"                     ~ 2,
      freq_sessoes == "Mensal"                         ~ 3,
      freq_sessoes == "Quinzenal"                      ~ 4,
      freq_sessoes == "Semanal"                        ~ 5
    ),
    papel_mesa = relevel(factor(papel_mesa), ref = "Jogador")
  )

stopifnot(!is.ordered(df$compra_pdf))
stopifnot(sum(is.na(df[, c("renda_num", "freq_sessoes_num")])) == 0)

# 3.3 Modelo completo
modelo_completo <- glm(
  consumidor ~ escala_valorizacao_posse + escala_resistencia_digital +
    renda_num + freq_sessoes_num + papel_mesa + compra_pdf,
  data   = df,
  family = binomial(link = "logit")
)
summary(modelo_completo)
vif(modelo_completo)

# 3.4 Seleção stepwise e modelo final
modelo_step <- stepAIC(modelo_completo, direction = "both", trace = TRUE)
summary(modelo_step)

anova(modelo_step, modelo_completo, test = "Chisq")

# 3.5 Qualidade de ajuste
hoslem.test(modelo_step$y, fitted(modelo_step), g = 10)

r2(modelo_step)
1 - (modelo_step$deviance / modelo_step$null.deviance)

# 3.6 Classificação (ponto de corte 0,50)
prob_predita   <- predict(modelo_step, type = "response")
classe_predita <- if_else(prob_predita >= 0.5, 1, 0)
matriz_confusao <- table(Observado = df$consumidor, Predito = classe_predita)
print(matriz_confusao)

acuracia       <- sum(diag(matriz_confusao)) / sum(matriz_confusao)
sensibilidade  <- matriz_confusao[2, 2] / sum(matriz_confusao[2, ])
especificidade <- matriz_confusao[1, 1] / sum(matriz_confusao[1, ])
cat("Acurácia:", round(acuracia, 3),
    " | Sensibilidade:", round(sensibilidade, 3),
    " | Especificidade:", round(especificidade, 3), "\n")

# 3.7 Curva ROC, AUC e coeficiente de Gini
roc_obj <- roc(df$consumidor, prob_predita,
               levels = c("Não comprador", "Comprador"), direction = "<")
plot(roc_obj, print.auc = TRUE)

auc_in  <- as.numeric(auc(roc_obj))
gini_in <- 2 * auc_in - 1
round(c(AUC = auc_in, Gini = gini_in), 4)

# 3.8 Validação cruzada (k-fold)
df <- df %>%
  mutate(consumidor_cv = factor(consumidor,
                                levels = c("Não comprador", "Comprador"),
                                labels = c("Nao_Comprador", "Comprador")))

set.seed(123)
ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE,
                     summaryFunction = twoClassSummary)
modelo_cv <- train(
  update(formula(modelo_step), consumidor_cv ~ .),
  data = df, method = "glm", family = "binomial",
  trControl = ctrl, metric = "ROC"
)
print(modelo_cv)

auc_cv  <- modelo_cv$results$ROC
gini_cv <- 2 * auc_cv - 1
round(c(AUC_cv = auc_cv, Gini_cv = gini_cv), 4)

# 3.9 Odds ratios e intervalos de confiança
odds_ratios <- exp(cbind(OR = coef(modelo_step), confint(modelo_step)))
print(round(odds_ratios, 3))


###############################################################################
# ETAPA 4 — Análises complementares
###############################################################################

# 4.1 Sensibilidade: modelo sem o construto de resistência à substituição digital
modelo_sem_f2 <- glm(
  consumidor ~ escala_valorizacao_posse + renda_num + freq_sessoes_num +
    papel_mesa + compra_pdf,
  data   = df,
  family = binomial(link = "logit")
)
summary(modelo_sem_f2)

modelo_sem_f2_step <- stepAIC(modelo_sem_f2, direction = "both", trace = FALSE)
summary(modelo_sem_f2_step)

roc_sem_f2 <- roc(df$consumidor, predict(modelo_sem_f2_step, type = "response"),
                  levels = c("Não comprador", "Comprador"), direction = "<")
round(as.numeric(auc(roc_sem_f2)), 3)

# 4.2 Taxa de referência e pontos de corte alternativos
taxa_referencia <- mean(df$consumidor == "Comprador")
round(100 * taxa_referencia, 1)

youden <- roc_obj$sensitivities + roc_obj$specificities - 1
corte_youden      <- roc_obj$thresholds[which.max(youden)]
corte_prevalencia <- sort(prob_predita, decreasing = TRUE)[round(taxa_referencia * nrow(df))]

metricas_corte <- function(corte) {
  previsto  <- prob_predita >= corte
  observado <- df$consumidor == "Comprador"
  c(corte          = corte,
    acuracia       = mean(previsto == observado),
    sensibilidade  = sum(previsto & observado) / sum(observado),
    especificidade = sum(!previsto & !observado) / sum(!observado))
}
round(rbind(corte_050    = metricas_corte(0.50),
            youden       = metricas_corte(corte_youden),
            prevalencia  = metricas_corte(corte_prevalencia)), 3)

# 4.3 Validação cruzada aninhada (seleção stepwise refeita em cada partição)
set.seed(123)
n_repeticoes <- 20

auc_teste <- function(modelo, teste) {
  as.numeric(auc(roc(teste$consumidor, predict(modelo, newdata = teste, type = "response"),
                     levels = c("Não comprador", "Comprador"), direction = "<")))
}

resultado_cv <- data.frame()
selecoes <- character()
for (r in seq_len(n_repeticoes)) {
  folds <- createFolds(df$consumidor, k = 10)
  for (idx in folds) {
    treino <- df[-idx, ]
    teste  <- df[idx, ]
    m_fixo     <- glm(formula(modelo_step), data = treino, family = binomial)
    m_aninhado <- stepAIC(glm(formula(modelo_completo), data = treino, family = binomial),
                          direction = "both", trace = FALSE)
    resultado_cv <- rbind(resultado_cv,
                          data.frame(auc_formula_fixa = auc_teste(m_fixo, teste),
                                     auc_aninhada     = auc_teste(m_aninhado, teste)))
    selecoes <- c(selecoes,
                  paste(sort(attr(terms(m_aninhado), "term.labels")), collapse = " + "))
  }
}
round(colMeans(resultado_cv), 3)

termos <- c("escala_valorizacao_posse", "escala_resistencia_digital", "renda_num",
            "freq_sessoes_num", "papel_mesa", "compra_pdf")
round(sapply(termos, function(t) mean(grepl(t, selecoes, fixed = TRUE))), 2)
round(sort(table(selecoes) / length(selecoes), decreasing = TRUE), 2)

# 4.4 Robustez da codificação e covariáveis adicionais
df <- df %>%
  mutate(renda_f        = factor(renda, ordered = FALSE),
         genero_masc    = as.numeric(genero == "Masculino"),
         tempo_rpg_num  = as.numeric(tempo_rpg))

modelo_renda_dummy <- update(modelo_completo, . ~ . - renda_num + renda_f)
anova(modelo_completo, modelo_renda_dummy, test = "Chisq")
AIC(modelo_completo, modelo_renda_dummy)

for (v in c("idade", "genero_masc", "tempo_rpg_num")) {
  m <- update(modelo_step, as.formula(paste(". ~ . +", v)))
  cat("\n+", v, "| AIC:", round(AIC(m), 2), "(modelo final:", round(AIC(modelo_step), 2), ")\n")
  print(round(summary(m)$coefficients[v, ], 4))
}

# 4.5 Calibração do modelo final e forma funcional
hl <- hoslem.test(modelo_step$y, fitted(modelo_step), g = 10)
round(cbind(hl$observed, hl$expected), 2)

for (g in c(6, 8, 10, 12, 15)) {
  h <- hoslem.test(modelo_step$y, fitted(modelo_step), g = g)
  cat("g =", g, "| chi2 =", round(h$statistic, 2), "| p =", round(h$p.value, 4), "\n")
}

df <- df %>%
  mutate(val_c = escala_valorizacao_posse - mean(escala_valorizacao_posse),
         res_c = escala_resistencia_digital - mean(escala_resistencia_digital))

modelo_forma <- update(modelo_step, . ~ . + I(val_c^2) + I(res_c^2) + val_c:res_c)
anova(modelo_step, modelo_forma, test = "Chisq")

