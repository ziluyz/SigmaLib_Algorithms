# SigmaLib_Algorithms (crate `algorithms`)

Библиотека переиспользуемых численных алгоритмов для вычислительных проектов Sigma: собственные
реализации на Rust и FFI-мост к GNU Scientific Library. От `processor` и Sigma не зависит.

Общая картина — `../CLAUDE.md`. План тестов — раздел «Система тестов» в `../ROADMAP.md`.

Имя папки и имя crate различаются: в зависимых проектах
`algorithms = { path = "../SigmaLib_Algorithms" }`.

## Команды

```bash
cargo build
cargo test
cargo test -- --ignored     # тесты известных дефектов (сейчас падают — так и задумано)
```

## Тесты

Интеграционные тесты в `tests/`: `algebra.rs`, `argument_root.rs`, `gsl_bessel.rs`,
`gsl_integration.rs`. Сравнение чисел — только макросом `assert_close!(a, b, abs = …, rel = …)`
из `tests/common/mod.rs` (подключается строкой `mod common;`), допуск указывается явно.
Эталонные значения сопровождать источником. Тест с `#[ignore = "known bug: …"]` описывает
правильное поведение, которого пока нет: исправление = снять `#[ignore]`.

`find_zeros_by_argument` засоряет вывод тестов отладочными `println!` — смотреть результат через
`cargo test 2>&1 | grep "test result"`.

Требуется системный GSL: `sudo apt install libgsl-dev`. `build.rs` линкует `gsl` и `gslcblas`
напрямую (`cargo:rustc-link-lib`), без pkg-config.

## Состав

| Модуль | Содержимое |
|---|---|
| `algebra` | `find_root` (бисекция), `find_root_sequence` (n последовательных корней от стартовой точки, точность 1e-7) |
| `algebra::argument_root` | `find_zeros_by_argument` — нули аналитической функции в прямоугольнике комплексной плоскости по принципу аргумента с адаптивным дроблением (`RectF64`, `Limits`, `Root`, `StopReason`, `EvalCache`). По git-истории — «beginning», не завершён; содержит отладочные `println!` |
| `gsl::errors` | `turn_off_error_handler()`, `get_error_description(status)` |
| `gsl::integration` | `qng` (функция), `Qawo` (осциллирующие интегранды: `new` → `initialize` → `integrate`, `change_parameters`/`change_limits`), `Romberg` (`new` → `initialize` → `integrate`) |
| `gsl::sfunctions::bessel` | `jn(n, x)`, `djn(n, x)` (производная через рекуррентность) |

## Соглашения

- Все публичные функции возвращают `Result<_, String>`; функции-аргументы у `algebra` —
  `Fn(f64) -> Result<f64, String>`, у интеграторов GSL — `FnMut(f64) -> f64`.
- FFI: объявления `extern "C"` лежат в `gsl_signatures.rs` рядом с модулем-обёрткой; наружу
  торчат только безопасные обёртки. Замыкание передаётся в GSL через `GslFunction` +
  `trampoline` (`gsl/integration.rs`): в `params` кладётся указатель на `&mut dyn FnMut`.
- Ресурсы GSL (workspace, таблицы) оборачиваются в структуру с `Drop`.
- Обработчик ошибок GSL по умолчанию вызывает `abort()`. Коды возврата проверяются обёртками,
  но это работает только после `gsl::errors::turn_off_error_handler()` — вызывающий проект
  обязан вызвать его один раз в начале расчёта.
- Документация функций — doc-комментарии с разделами Arguments / Result / Errors (как в `algebra.rs`).
- Библиотека не печатает в stdout и не паникует на пользовательских данных.
- Новый алгоритм добавляется вместе с тестами: аналитический случай + граничные/ошибочные входы.

## Известные проблемы

- `Romberg::drop` пуст — workspace не освобождается; повторный `initialize` у `Romberg`/`Qawo` течёт.
- `bessel::jn` не проверяет статус GSL (используется не-`_e` вариант функции).
- `find_root_sequence` зацикливается навсегда, если корней меньше `n` (нет ограничения поиска)
  или если функция возвращает ошибку (она неотличима от «нет корня на отрезке»). Закреплено
  двумя помеченными тестами.
- `find_zeros_by_argument`: для `z² + 1` находит один корень из двух; на кратном корне
  завершается ошибкой «Contour too close to a zero». Закреплено двумя помеченными тестами.
- Паника внутри замыкания в `trampoline` раскручивается через кадры C — неопределённое поведение.
- `qawo`: ошибка округления GSL намеренно игнорируется (см. коммит `ac3c72e`).
