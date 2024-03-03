# DA TRADURRE:



#Funzioni aggiuntive da implementare:
# - psn
# - psc
# - pst

function (x, xi = 0, omega = 1, alpha = 0, nu = Inf, dp = NULL, 
          method = 0, lower.tail = TRUE, log.p = FALSE, ...) 
{
  if (!is.null(dp)) {
    if (!missing(alpha)) 
      stop("You cannot set both component parameters and dp")
    xi <- dp[1]
    omega <- dp[2]
    alpha <- dp[3]
    nu <- dp[4]
  }
  if (length(alpha) > 1) 
    stop("'alpha' must be a single value")
  if (length(nu) > 1) 
    stop("'nu' must be a single value")
  if (nu <= 0) 
    stop("'nu' must be positive")
  dp.std <- c(0, 1, alpha, nu)
  delta <- alpha/sqrt(1 + alpha^2)
  if (nu == Inf) 
    return(psn(x, xi, omega, alpha))
  if (nu == 1) 
    return(psc(x, xi, omega, alpha))
  int.nu <- (round(nu) == nu)
  if (method < 0 | method > 5 | method != round(method)) 
    stop("invalid 'method' value")
  if ((method == 1 | method == 4) & !int.nu) 
    stop("selected 'method' does not work for non-integer nu")
  z <- (x - xi)/omega
  pr <- rep(NA, length(z))
  ok <- !(is.na(z) | (z == Inf) | (z == -Inf) | (omega <= 0))
  z <- z[ok]
  nu0 <- (8.2 + 3.55 * log(log(length(z) + 1)))
  if (alpha == 0) 
    p <- pt(z, df = nu)
  else if (abs(alpha) == Inf) {
    z0 <- replace(z, alpha * z < 0, 0)
    p <- pf(z0^2, 1, nu)
    if (alpha < 0) 
      p <- (1 - p)
  }
  else {
    fp <- function(v, alpha, nu, t.value) psn(sqrt(v) * t.value, 
                                              0, 1, alpha) * dchisq(v * nu, nu) * nu
    if (method == 4 || (method == 0 && int.nu && (nu <= nu0))) {
      p. <- pst_int(z, 0, 1, alpha, nu)
      p <- if (lower.tail) 
        p.
      else 1 - p.
      p <- if (log.p) 
        log(p)
      else p
    }
    else {
      p <- numeric(length(z))
      for (i in seq_len(length(z))) {
        if (abs(z[i]) == Inf) 
          p[i] <- (1 + sign(z[i]))/2
        if (method == 5 | method == 0 & abs(z[i]) > (30 + 
                                                     1/sqrt(nu))) {
          lp <- st_tails(z[i], alpha, nu, lower.tail = lower.tail)
          p[i] <- if (log.p) {
            if (z[i] < 0) 
              lp
            else log(1 - exp(lp))
          }
          else {
            if (z[i] < 0) 
              exp(lp)
            else 1 - exp(lp)
          }
        }
        else {
          if (method == 1 || (method == 0 && int.nu && 
                              (nu > nu0))) {
            out <- try(pmst(z[i], 0, matrix(1, 1, 1), 
                            alpha, nu, ...), silent = TRUE)
            p. <- if (inherits(out, "try-error")) 
              NA
            else out
          }
          else {
            upper <- 10 + 50/nu
            if (method == 2 || (method == 0 & (z[i] < 
                                               upper))) {
              p0 <- acos(delta)/pi
              int <- integrate(dst, min(0, z[i]), max(0, 
                                                      z[i]), dp = dp.std, stop.on.error = FALSE, 
                               ...)
              p. <- p0 + sign(z[i]) * int$value
            }
            else {
              p. <- integrate(fp, 0, Inf, alpha, nu, 
                              z[i], stop.on.error = FALSE, ...)$value
            }
          }
          p[i] <- if (lower.tail) 
            p.
          else 1 - p.
          p[i] <- if (log.p) 
            log(p[i])
          else max(0, min(1, p[i]))
        }
      }
    }
  }
  pr[ok] <- p
  pr[x == Inf] <- if (log.p) 
    0
  else 1
  pr[x == -Inf] <- if (log.p) 
    -Inf
  else 0
  pr[omega <= 0] <- NaN
  names(pr) <- names(x)
  return(pr)
}
<bytecode: 0x00000200b328e1d0>
  <environment: namespace:sn>